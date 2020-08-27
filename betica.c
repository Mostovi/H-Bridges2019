#include "Header.h"

//#include <string.h>
//Measurements

/********************************** definisanje MAKROa koji se koriste u proracunima */
#define SAMPLING_TIME 0.546134
#define SAMPLING_FREQ (480000/16)
#define NUM_OF_POINTS (SAMPLING_FREQ*SAMPLING_TIME)
#define F_FUNDAMENTAL 10 //ZA POCETNO POGADJANJE UCESTANOSTI RSH
#define F0_FILTER 230   //ZA POCETNO POGADJANJE UCESTANOSTI RSH
#define DELTA_W_BPF 85   //3dB PROPUSNI OPSEG BANDPASS FILTRA
#define KSI_PLL 9		 //KOEFICIJENT PRIGUSENJA PLL PETLJE
#define WN_PLL 150		 //PRIRODNA UCESTANOST PLL PETLJE
#define R_NOTCH 0.8      //DEFINISE POLOVE I NULE NOTCH FILTRA. ZA 1 JE TEORIJSKI BESKONACNO SLABLJENJE ALI JE SPORIJI FILTAR. OVO 0.8 JE VEOMA BRZO A SLABLJENJE NE TREBA DA BUDE VECE OD TOGA JER I SAM PLL RADI KAO LPF
#define KP_OBS 0.01		 //PROPORCIONALNO POJACANJE OBSERVERA ZA TAU_OBS = 12.5 ms
#define KI_OBS 15		 //INTEGRALNO POJACANJE OBSERVERA ZA TAU_OBS = 67 ms
#define KI_OBS_BPF 80    //INTEGRALNO POJACANJE OBSERVERA ZA TAU_OBS = 12.5 ms (OVAJ OBSERVER SE KORISTI ZA ADAPTIRANJE CENTRALNE UCESTANOSTI BANDPASS I NOTCH FILTARA)
#define SCALE_FACTOR 60  //PRETPOSTAVLJENA (IZ MATLABA) POCETNA AMPLITUDA IZLAZNOG SIGNALA IZ BANDPASS FILTRA (ZLEBNOG HARMONIKA). OVA VREDNOST SE POSLE MENJA IZLAZOM IZ DEMODULATORA

#pragma DATA_SECTION(Ia, "DMARAML4");
#pragma DATA_SECTION(Einverse, "DMARAML4");
#pragma DATA_SECTION(Ualpha, "DMARAML4");
#pragma DATA_SECTION(Ubeta, "DMARAML4");
#pragma DATA_SECTION(_sin, "DMARAML4");
#pragma DATA_SECTION(_cos, "DMARAML4");
#pragma DATA_SECTION(Ud, "DMARAML4");
#pragma DATA_SECTION(Uq, "DMARAML4");
#pragma DATA_SECTION(Ua, "DMARAML4");
#pragma DATA_SECTION(Ub, "DMARAML4");
#pragma DATA_SECTION(Uc, "DMARAML4");
#pragma DATA_SECTION(da, "DMARAML4");
#pragma DATA_SECTION(db, "DMARAML4");
#pragma DATA_SECTION(dc, "DMARAML4");
#pragma DATA_SECTION(DMABuff, "DMARAML4");
#pragma DATA_SECTION(Measurements, "DMARAML4");



Uint16 ErrorFlagTZ= 0;
Uint16 dmaCnt = 0;
volatile Uint16 DMABuff[16];
Uint32 Measurements[6] = {0, 0, 0, 0, 0, 0};
float32 Ia = 0.0f;
float32 omegaFiltered = 0, omegaTestMax2 = 10.0f,omegaTestMax = 10.0f, omegaTest = 6.28f, thetaTest = 0.0f,omegaIncrement=0.0001f;
float32 _sin[3] = {0.0f, 0.0f};
float32 _cos[3] = {0.0f, 0.0f};
//Voltages
float32 Ualpha[2] = {0.0f,0.0f}, Ubeta[2] = {0.0f,0.0f}, Ua = 0.0f, Ub = 0.0f, Uc = 0.0f;
float32 kUd = 2.0f, Ud = 1.0f, Uq = 0.0f;
float32 E = 190.0f, Einverse = 0.005263f;
//SVPWM
Uint16 da = 0, db =0, dc = 0;
float32 dta = 0.0f, dtb = 0.0f, dtc = 0.0f;
//Flags
Uint16 RegulationEnabled = 0, EstimationEnabled=0;
//Data Output
Uint16 outerDataCount, dataCount = 9999, outerDataMaxCount = 100;
float32 outData[1536] = {};

float mains_frequency;
long int j=0,dma_counter=0,brojac=0, cos_table_index=0,speed_pointer=0,cos_table_index_sin=7500,speed_pointer_lag,counter_remaining=12400;
int k=0,i1,kk,transient=0,speed_diff;
volatile float y4;
volatile float32 I_line;
volatile float32 brzina;

//promenljive za bpf
int f0_filter,N_kask;
float32 delta_w_kask,delta_w_stop_kask,ks,qs,q,delta_w_pojedinacnog,delta_w_rel,w0_rel,alfa,beta,one_minus_alfa_beta,two_alfa_minus,y1,y1_bef,y1_befbef,sum_befbef1,sum_bef1,
y2,y2_bef,y2_befbef,sum_befbef2,sum_bef2,y3,y3_bef,y3_befbef,sum_befbef3,sum_bef3,y4_bef,y4_befbef,sum_befbef4,sum_bef4;

//promenljive za PLL
volatile float32 theta_pll,e, e_filt,e_filt_preth, e_filt_bef, e_filt_befbef, e_bef, e_befbef, w_0_notch_rel,w_0_notch_rel_const,b1_notch,r,a1_notch,a2_notch,
e_bef_demod,e_befbef_demod,e_filt_demod,e_filt_bef_demod,e_filt_befbef_demod,RSH_demod,error_demod_obs,error_obs_zajedno,RSH_demod_obs,ki_demod_zajedno,
 ki_pll, kp_pll,kp_pll_coeff ,ki_pll_coeff, Ts_pll, VCO, up_pll, ui_pll, ki_pll_zajedno, f_pll, theta_increment,two_pi,error_obs,f_obs,f_obs_bpf,up_obs,ui_obs,
 ui_obs_bpf,kp_obs,ki_obs,ki_obs_zajedno,error_obs_bpf,ki_obs_zajedno,ki_obs_bpf,ki_obs_bpf_zajedno,theta_table_incr=2*PI/10000,
fs=150000000/(156*2*16);// TODO JEBI SE TACNA UCESTANOST

POSSPEED qep_posspeed=POSSPEED_DEFAULTS;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void main(void)
{
    //BUFFER DISABLE
    GpioDataRegs.GPASET.bit.GPIO31 = 1;
    GpioDataRegs.GPASET.bit.GPIO23 = 1;     //5V HOLD

    DSP_INIT_CONFIG();
    forceOneShotPwm();

    EnablePeripheralInterrupts();
    //disableTrans();
    StartPeripherals();
    DELAY_US(20);
    ConfigBPF();
    ConfigPLL();
    sincos(0.0f, &_sin[0], &_cos[0]);
    qep_posspeed.init(&qep_posspeed);
    //BUFFER ENABLE

    qep_posspeed.mech_scaler = 1/4096.0f;
    qep_posspeed.pole_pairs = 2;
    GpioDataRegs.GPACLEAR.bit.GPIO31 = 1;
    int ii=0;
    int i;
    for(i=0;i<1000;i++)
    {
      ii++;
    }
    GpioDataRegs.GPASET.bit.GPIO31 = 1;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ResetRegulationVars(void)
{
    ConfigBPF();
    ConfigPLL();

    omegaTest = 0.0f;
    Ud = 0.0f;
}

//  ?  // END REGULATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__interrupt void
DmaInterrupt(void)
{
     GpioDataRegs.GPBSET.bit.GPIO60 = 1;
     ReadMeasurementBuffer();
     GpioDataRegs.GPBCLEAR.bit.GPIO60 = 1;
     PieCtrlRegs.PIEACK.all = PIEACK_GROUP7;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__interrupt void
InverterInterrupt(void)
{
    GpioDataRegs.GPBSET.bit.GPIO62 = 1;
    //I_line = (float32)(Measurements[0]);
    I_line = (float32)(Measurements[0]>>1);
   // I_line = (I_line - 2509.0f)*1.0f; // TODO ZAMENITI 1.0F SA KOEFICIJENTOM SKALIRANJA REALNE STRUJE
    Measurements[0] = 0;
    
    if (RegulationEnabled)
    {
        // ZADAVANJE NAPONA
        if(omegaTest<omegaTestMax)omegaTest += omegaIncrement;
        else if(omegaTest>omegaTestMax)omegaTest -= omegaIncrement;
        if(omegaTest<1.0f)omegaTest = 1.0f;

        thetaTest+=omegaTest*DVAPI*0.000033f;
        if(thetaTest>DVAPI){thetaTest-=DVAPI;}
        if(thetaTest<0.0f){thetaTest+=DVAPI;}

        Ud = kUd * omegaTest ; // V/F
        Uq = 0;

        sincos(thetaTest, &_sin[0], &_cos[0]);
        Ualpha[0] = Ud * ( _cos[0]  ) - Uq * ( _sin[0] ) ;
        Ubeta[0]  = Ud * ( _sin[0]  ) + Uq * ( _cos[0] ) ;

        // * Inverse Clarke's transformation
        Ua = Ualpha[0];
        Ub = (1.73205081f * Ubeta[0] - Ualpha[0]) * 0.5f;
        Uc = -(1.73205081f * Ubeta[0] + Ualpha[0]) * 0.5f;
        // modulation signals
        da = (Uint16)(INVERTER_PRD * (0.5f - Ua * Einverse) );//- dta); // minus, ne dirati!!!
        db = (Uint16)(INVERTER_PRD * (0.5f - Ub * Einverse) );//- dtb);
        dc = (Uint16)(INVERTER_PRD * (0.5f - Uc * Einverse) );//- dtc);
        // Modulation signals limiter
        if(da>INV_MAX_D)da=INV_MAX_D;else if (da<INV_MIN_D) da=INV_MIN_D;
        if(db>INV_MAX_D)db=INV_MAX_D;else if (db<INV_MIN_D)db=INV_MIN_D;
        if(dc>INV_MAX_D)dc=INV_MAX_D;else if (dc<INV_MIN_D)dc=INV_MIN_D;

        PwmInverter_1.CMPA.half.CMPA = da + DEADTIMEINVERTER_POLA;
        PwmInverter_1.CMPB = da - DEADTIMEINVERTER_POLA;
        PwmInverter_2.CMPA.half.CMPA = db + DEADTIMEINVERTER_POLA;
        PwmInverter_2.CMPB = db - DEADTIMEINVERTER_POLA;
        PwmInverter_3.CMPA.half.CMPA = dc + DEADTIMEINVERTER_POLA;
        PwmInverter_3.CMPB = dc - DEADTIMEINVERTER_POLA;
        GpioDataRegs.GPBCLEAR.bit.GPIO62 = 1;
    }
    else {
        dataCount =0;
        ResetRegulationVars();
        disableTrans();
    }

    GpioDataRegs.GPBSET.bit.GPIO62 = 1;
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    
    /* Bandpass filtriranje */
	/*******************************************************************************************************************************************************************************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
	// SLUZI DA IZ CELOG SPEKTRA STATORSKE STRUJE (PROMENLJIVA I_line) IZDVOJI SIGNAL U OKOLINI SVOJE CENTRALNE UCESTANOSTI (KOJA SE ADAPTIRA TAKO DA TO UVEK BUDE ZLEBNI HARMONIK)
	// FILTAR JE REALIZOVAN KASKADNO PRI CEMU SVAKA DODATNA KASKADA JOS VISE FILTRIRA SIGNAL (STRMIJA JE KARAKTERISTIKA EKVIVALENTNOG FILTRA). PARAMETRI SVE CETIRI KASKADE SU ISTI I TE PROMENLJIVE SU GRUPISANE U OKVIRU FUNKCIJE void ConfigBPF(void) DA BI SE KORISTIO STO JE MANJE MOGUCI BROJ RACUNSKIH OPERACIJA
	// PARAMETAR BETA JE VEZAN ZA CENTRALNU UCESTANOST FILTRA A ALFA ZA SIRINU PROPUSNOG OPSEGA (BETA SE MENJA U TOKU RADA A ALFA NE)
	// U OVOM DELU KODA SVE PROMENLJIVE SU SKALARI. ZBOG DIFERENCNE PRIRODE JEDNACINA POTREBNO JE CUVATI ODREDJENI BROJ PRETHODNIH VREDNOSTI
	// PROMENLJIVE OBLIKA y SU IZLAZI IZ FILTARSKIH KASKADA U TRENUTKU [n]
	// PROMENLJIVE OBLIKA y_bef SU VREDNOSTI IZLAZA IZ PRETHODNE ITERACIJE [n-1]
	// PROMENLJIVE OBLIKA y_befbef SU VREDNOSTI IZLAZA OD PRE 2 ITERACIJE [n-2]
	// PROMENLJIVE OBLIKA sum_bef SU ULAZI U FILTARSKE KASKADE IZ PRETHODNE ITERACIJE [n-1]
	// PROMENLJIVE OBLIKA sum_befbef SU ULAZI U FILTARSKE KASKADE OD PRE 2 ITERACIJE [n-2]
	// IZLAZ SVAKE KASKADE JE ULAZ U NAREDNU (OSIM POSLEDNJE, CIJI JE IZLAZ ULAZ U PLL I PREDSTAVLJA ISFILTRIRANI SIGNAL RSH- y4)

    if(EstimationEnabled)
    {
        // Prva kaskada
        y1=one_minus_alfa_beta*y1_bef+two_alfa_minus*y1_befbef+alfa*(I_line-sum_befbef1);
        y1_befbef=y1_bef;
        sum_befbef1=sum_bef1;
        y1_bef=y1;
        sum_bef1=I_line;


        // Druga kaskada
        y2=one_minus_alfa_beta*y2_bef+two_alfa_minus*y2_befbef+alfa*(y1-sum_befbef2);
        y2_befbef=y2_bef;
        sum_befbef2=sum_bef2;
        y2_bef=y2;
        sum_bef2=y1;

        // Treca kaskada
        y3=one_minus_alfa_beta*y3_bef+two_alfa_minus*y3_befbef+alfa*(y2-sum_befbef3);
        y3_befbef=y3_bef;
        sum_befbef3=sum_bef3;
        y3_bef=y3;
        sum_bef3=y2;

        // Cetvrta kaskada
        y4=one_minus_alfa_beta*y4_bef+two_alfa_minus*y4_befbef+alfa*(y3-sum_befbef4);
        y4_befbef=y4_bef;
        sum_befbef4=sum_bef4;
        y4_bef=y4;
        sum_bef4=y3;
        //KRAJ BANDPASS FILTRIRANJA
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/


        /* Phase-locked loop algoritam */
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        // ULAZ JE IZLAZ POSLEDNJE KASKADE BANDPASS FILTRA
        // DETEKTOR FAZE NA BAZI MNOZACA, ULAZNA SINUSOIDA SE MNOZI IZLAZOM PLL-A KOJI U USTALJENOM STANJU TREBA DA BUDE KOSINUSOIDA (PI/2 FAZNI STAV U ODNOSU NA ULAZNI SIGNAL)
        // cos_table_index JE INDEX U NIZU KOSINUSNE TABELE PREDEFINISANE NA POCETKU PROGRAMA. TAJ INDEX SE POMERA NA KRAJU OVE INTERAPT RUTINE KADA SE PRORACUNA NOVI FAZNI STAV ULAZNOG SIGNALA (TAJ UGAO JE IZLAZ IZ PLL- A)
        // IZLAZ DETEKTORA FAZE JE PROMENLJIVA e, KOJA SADRZI I SIGNAL NA DVOSTRUKOJ UCESTANOSTI
        // TAJ SIGNAL DVOSTRUKE UCESTANOSTI SE FILTRIRA POMOCU NOTCH FILTRA
        // PROMENLJIVE U OKVIRU NOTCH FILTRA SU:
        // e_filt - IZLAZ IZ NOTCH FILTRA. TAJ SIGNAL SE KORISTI U DALJEM PROCESIRANJU SIGNALA KAO SIGNAL FAZNE RAZLIKE ZLEBNOG HARMONIKA I IZLAZA IZ PLL- A
        // e_filt_bef - IZLAZ IZ NOTCH FILTRA U TRENUTKU [n-1]
        // e_filt_befbef - IZLAZ IZ NOTCH FILTRA U TRENUTKU [n-2]
        // e - IZLAZ DETEKTORA FAZE, ULAZ U PLL U TRENUTKU [n]. PORED FAZNE RAZLIKE OVAJ SIGNAL IMA I POJACANJE JEDNAKO POLOVINI AMPLITUDE SIGNALA KOJI PREDSTAVLJA ZLEBNI HARMONIK. TO SE KOMPENZUJE PROMENOM POJACANJA kp I ki POMOCU DEMODULATORA (KASNIJE)
        // e_bef - ULAZ U PLL U TRENUTKU [n-1]
        // e_befbef - ULAZ U PLL U TRENUTKU [n-2]
        // NAKON DETEKTORA FAZE, NJEGOV IZLAZ SE SALJE U PI REGULATOR KOJI IMA CILJ DA GRESKU IZMEDJU FAZNOG STAVA DVA SIGNALA POSTAVI NA 0 (ODNOSNO PI/2 AKO SE SMATRA DA SU OBE SINUSOIDE A NE SIN I COS)

        //DETEKTOR FAZE
        // lol
        sincos(theta_pll, &_sin[0], &_cos[0]);
        e=_cos[0]*y4;	// MNOZENJE DVE SINUSOIDE

        e_filt=-a1_notch*e_filt_bef-a2_notch*e_filt_befbef+e+b1_notch*e_bef+e_befbef; // NOTCH FILTAR ZA DVOSTRUKU UCESTANOST
        e_befbef=e_bef;
        e_bef=e;
        e_filt_befbef=e_filt_bef;
        e_filt_bef=e_filt;

        //PI REGULATOR, ULAZ JE FILTRIRANI IZLAZ FAZNOG DETEKTORA (THETA_RSH-THETA_PLL)*0.5*I_RSH_MAX.
        //UTICAJ AMPLITUDE RSH SE ELIMINISE POMOCU DEMODULATORA (MENJAJU SE POJACANJA PI REGULATORA, MOGAO JE I ULAZNI SIGNAL DA SE DELI TIME, SVEJEDNO JE)
        up_pll=kp_pll*e_filt;
        ui_pll=ui_pll+ki_pll_zajedno*(e_filt+e_filt_preth);
        e_filt_preth=e_filt;

        //IZLAZ PI REGULATORA JE PRVO UCESTANOST RSH, KOJA SE ZATIM MNOZI SA 2PI I DODAJE U AKUMULATOR KOJI PREDSTAVLJA TRENUTNI UGAO
        //PROMENLJIVA ui_pll SADRZI POCETNU VREDNOST FO_FILTER RADI BRZEG USPOSTAVLJANJA STACIONARNOG STANJA. TO JE MOGLO DA SE ODRADI I OVDE (KAO ZA OBSERVER)
        f_pll=up_pll+ui_pll;

        //FIZICKA OGRANICENJA
        if (f_pll<0)
            f_pll=0;

        //GORNJI LIMIT POMAZE DA SE BRZE UDJE U STACIONARNO STANJE, ZAVISI OD RADNE UCESTANOSTI
        if (f_pll>1200)
            f_pll=1200;

        //INKREMENT UGLA KOJI PREDSTAVLJA IZLAZ PLL-A
        theta_pll+=theta_increment*f_pll;

        //OGRANICAVANJE PROMENLJIVE UGLA NA 2*PI
        if (theta_pll>=two_pi)
        {
            theta_pll-=two_pi;
        }

        //DEMODULACIJA AMPLITUDE RSH
        //IZLAZ DEMODULATORA JE 0.5*I_RSH_MAX.
        RSH_demod=_sin[0]*y4; //OVAJ SIGNAL TAKODJE SADRZI KOMPONENTU DVOSTRUKE UCESTANOSTI KOJU FILTRIRAMO ISTIM NOTCH FILTROM KAO I FAZNI STAV

        //NOTCH FILTAR
        //PROMENLJIVA e_filt_demod PREDSTAVLJA 0.5*I_RSH_MAX BEZ BRZOPROMENLJIVE KOMPONENTE. UPRAVO JE TO VELICINA KOJOM TREBA PODELITI POJACANJA kp I ki DA BI PETLJA RADILA SA KONSTANTNIM POJACANJIMA NEZAVISNO OD AMPLITUDE ULAZNOG SIGNALA KOJI TREBA PRATITI.
        e_filt_demod=-a1_notch*e_filt_bef_demod-a2_notch*e_filt_befbef_demod+RSH_demod+b1_notch*e_bef_demod+e_befbef_demod;
        e_befbef_demod=e_bef_demod;
        e_bef_demod=RSH_demod;
        e_filt_befbef_demod=e_filt_bef_demod;
        e_filt_bef_demod=e_filt_demod;

        //POSTAVLJA SE DONJI LIMIT DA BI SE OBEZBEDILO OD NUMERICKE GRESKE UKOLIKO SE U NEKOJ PREKIDNOJ RUTINI JAVI MNOGO MALI IZLAZ DEMODULATORA
        if (e_filt_demod>30)
        {
            kp_pll=kp_pll_coeff/e_filt_demod;
            ki_pll_zajedno=(ki_pll_coeff/e_filt_demod)*0.5*Ts_pll;
        }
        //KRAJ PLL-A
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/


        /* OBSERVER */
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        //IMA DVE ULOGE
        //PRVA JE DA FILTRIRA IZLAZ IZ PLL- A DA ESTIMIRANA BRZINA NE BI SADRZILA PREVELIKI SUM
        //DRUGA JE DA ADAPTIRA BANDPASS I NOTCH FILTRE JER SE ONI NE MOGU DIREKTNO ADAPTIRATI IZ PLL-A JER IMAJU BRZU DINAMIKU
        //FORMA OBSERVERA JE CLOSED LOOP LPF U VIDU PI REGULATORA
        //VREMENSKA KONSTANTA FILTRA SE BIRA TAKO DA BUDE VECA OD TRAJANJA PRELAZNOG REZIMA BANDPASS FILTRA (DA SE NE BI JAVILA NESTABILNOST) A KRACA OD MEHANICKE VREMENSKE KONSTANTE MASINE (DA SE NE BI JAVILO KASNJENJE)

        //IMPLEMENTIRANA SU DVA OBSERVERA RADI POREDJENJA NJIHOVIH IZLAZA (JEDAN BRZI I JEDAN SPORIJI)

        //OVO JE SPORIJI BPF KOJI REZULTUJE VECIM KASNJENJEM ALI ODZIVOM SA VEOMA MALO SUMA
        error_obs=f_pll-f_obs;
        up_obs=kp_obs*error_obs;
        ui_obs+=ki_obs_zajedno*error_obs;
        f_obs=F0_FILTER+up_obs+ui_obs;

        //OVO JE BRZI BPF POMOCU KOJEG SE ADAPTIRA BPF DA ZLEBNI HARMONIK NE BI ISPAO IZ OPSEGA BPF USLED KASNJENJA ODZIVA OBSERVERA
        error_obs_bpf=f_pll-f_obs_bpf;
        ui_obs_bpf+=ki_obs_bpf_zajedno*error_obs_bpf;
        f_obs_bpf=F0_FILTER+kp_obs*error_obs_bpf+ui_obs_bpf;

        /* I NOTCH I BPF SE DALJE ADAPTIRAJU POMOCU IZLAZA DRUGOG OBSERVERA */

        //ADAPTIRANJE PARAMETARA NOTCH FILTRA
        w_0_notch_rel=w_0_notch_rel_const*f_obs_bpf;
        // lol w-nothc
        sincos(w_0_notch_rel, &_sin[1], &_cos[1]);
        b1_notch=-2*_cos[1];
        a1_notch=b1_notch*r;

        sincos(w_0_notch_rel/2, &_sin[2], &_cos[2]);
        //ADAPTIRANJE PARAMETARA BANDPASS FILTRA. CENTRALNA UCESTANOST JE JEDNAKA POLOVINI ONE OD NOTCH FILTRA.
        one_minus_alfa_beta=(1-alfa)*2*_cos[2];
        //KRAJ OBSERVER-A
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/


        /* PUNJENJE IZLAZNOG NIZA ZA PRIKAZ ESTIMACIJE */
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/

            //TODO ODJE DRUGA VERZIJA DIREKTNO PRIKAZUJE UCESTANOST RSH
            //speed_buffer[speed_pointer]=(Uint16)f_obs;
            //speed_buffer_unfiltered[speed_pointer]=(Uint16)f_pll;
            //speed_buffer_bpf[speed_pointer]=(Uint16)f_obs_bpf;
            //speed_pointer++;
    
        //KRAJ ESTIMATORA-A
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************************************************************************************************************************************************************************************************************************************************************************/
        printData();

    }
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    qep_posspeed.calc(&qep_posspeed);
    GpioDataRegs.GPBCLEAR.bit.GPIO62 = 1;
    InverterTimekeeper.ETCLR.bit.INT = 1;    // Clear epwm interrupt flag
    PieCtrlRegs.PIEACK.all |= PIEACK_GROUP3; // Acknowledge interrupt group to PIE
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printData()     //f. ispisivanje snimljenih rezultata
{
    if(dataCount==50)
    {
        omegaTestMax = omegaTestMax2;
    }
    outerDataCount++;
    if (outerDataCount >= outerDataMaxCount)
    {


        //outData[dataCount] = 1.363636f*(f_obs-F_FUNDAMENTAL);
        //outData[500 + dataCount] = 1.363636f*(f_pll-F_FUNDAMENTAL);
        //outData[1000 + dataCount] = omegaFiltered;

        omegaFiltered = 1.363636f*(f_obs_bpf-F_FUNDAMENTAL);

        outData[dataCount] = omegaFiltered;
        outData[500 + dataCount] = omegaTest;
        outData[1000 + dataCount] = qep_posspeed.SpeedRpm_fr;

        dataCount += dataCount > 499 ? 0 : 1;
        outerDataCount = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__interrupt void
AdcInterrupt(void)
{
    AdcRegs.ADCTRL2.bit.RST_SEQ1 = 1;        // Reset SEQ1
    AdcRegs.ADCST.bit.INT_SEQ1_CLR = 1;
    PieCtrlRegs.PIEACK.all |= PIEACK_GROUP1; // Acknowledge interrupt to PIE
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__interrupt void
TZINT(void)
{
    GpioDataRegs.GPBSET.bit.GPIO34 = 1;      //TURN OFF THE LED ON DSP CARD
    ErrorFlagTZ = 1;
    PieCtrlRegs.PIEACK.all |= PIEACK_GROUP2; // Acknowledge interrupt to PIE
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadMeasurementBuffer()     // ? // MEASUREMENTS
{
    for(dmaCnt=0;dmaCnt<16;dmaCnt++)
    {
        Measurements[0] += DMABuff[dmaCnt];    
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DisablePeripheralInterrupts()
{
    PieCtrlRegs.PIEIER3.bit.INTx1 = 0; // PWM 1 Interrupt
    PieCtrlRegs.PIEIER3.bit.INTx2 = 0; // PWM 1 Interrupt
    PieCtrlRegs.PIEIER3.bit.INTx4 = 0; // PWM 4 Interrupt

    PieCtrlRegs.PIEIER1.bit.INTx1 = 0; // SEQ1 Interrupt
    PieCtrlRegs.PIEIER1.bit.INTx6 = 0; // ADC Interrupt

    PieCtrlRegs.PIEACK.bit.ACK1 = 0;
    PieCtrlRegs.PIEIFR1.bit.INTx2 = 0;
    PieCtrlRegs.PIEIFR1.bit.INTx6 = 0;

    PieCtrlRegs.PIEIER7.bit.INTx1 = 0; // Dma Interrupt
    PieCtrlRegs.PIEIER2.bit.INTx1 = 0; // TZ EPWM1
    PieCtrlRegs.PIEIER1.bit.INTx4 = 0; // 1 // XINT1

    PieCtrlRegs.PIEIER1.bit.INTx7 = 0; //CPUTimer0
    PieCtrlRegs.PIEIER9.bit.INTx2 = 0; // 1 //SCIATX
    PieCtrlRegs.PIEIER9.bit.INTx1 = 0; // 1 //SCIARX
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EnablePeripheralInterrupts()
{
    PieCtrlRegs.PIEIER3.bit.INTx1 = 0;
    PieCtrlRegs.PIEIER3.bit.INTx2 = 0;
    PieCtrlRegs.PIEIER3.bit.INTx3 = 0;
    PieCtrlRegs.PIEIER3.bit.INTx4 = 1;
    PieCtrlRegs.PIEIER3.bit.INTx5 = 0;
    PieCtrlRegs.PIEIER3.bit.INTx6 = 0; // 1 // inverter timekeeper
    PieCtrlRegs.PIEIER1.bit.INTx1 = 0; // 1 SEQ1 Interrupt enable // bitno
    PieCtrlRegs.PIEIER1.bit.INTx6 = 0; // ADC Interrupt enable
    PieCtrlRegs.PIEIER1.bit.INTx8 = 0; // WAKEINT
    PieCtrlRegs.PIEACK.bit.ACK1 = 1;   // WAKEINT ACK
    PieCtrlRegs.PIEIER7.bit.INTx1 = 1; // 1 //Dma Interrupt enable
    PieCtrlRegs.PIEIER2.bit.INTx1 = 1; // 1 // TZ EPWM1
    PieCtrlRegs.PIEIER1.bit.INTx4 = 0; // XINT1

    //PieCtrlRegs.PIEIER1.bit.INTx7 = 1; //  1  //CPUTimer0
    PieCtrlRegs.PIEIER9.bit.INTx2 = 0; //  1  //SCIATX
    PieCtrlRegs.PIEIER9.bit.INTx1 = 0; //  1  //SCIARX
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void StartPeripherals()
{
    EALLOW;

    InverterTimekeeper.CMPA.half.CMPA = 3000; // dont interrupt
    BoostTimeKeeper.CMPA.half.CMPA = 2000;    // SYNC...
    BoostTimeKeeper.CMPB = 2000;

    ClearArray(DMABuff, 24);

    AdcRegs.ADCTRL2.bit.SOC_SEQ1 = 1;           //START RECEIVING SOC

    SysCtrlRegs.PCLKCR0.bit.ADCENCLK = 1;       //ENABLE PERIPHERAL CLOCKS
    SysCtrlRegs.PCLKCR0.bit.I2CAENCLK= 0;       //-
    SysCtrlRegs.PCLKCR0.bit.SCIAENCLK= 0;       //-
    SysCtrlRegs.PCLKCR0.bit.ECANAENCLK= 0;      //DISABLE PERIPHERAL CLOCKS  
    SysCtrlRegs.PCLKCR0.bit.ECANBENCLK= 0;      //- 
    SysCtrlRegs.PCLKCR0.bit.MCBSPAENCLK= 0;     //-   
    SysCtrlRegs.PCLKCR0.bit.MCBSPBENCLK= 0;     //-  
    SysCtrlRegs.PCLKCR0.bit.SCIBENCLK= 0;       //-
    SysCtrlRegs.PCLKCR0.bit.SCICENCLK= 0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP1ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP2ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP3ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP4ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP5ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.ECAP6ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.EPWM6ENCLK=0;       //-
    SysCtrlRegs.PCLKCR1.bit.EQEP1ENCLK=1;       //-
    SysCtrlRegs.PCLKCR1.bit.EQEP2ENCLK=0;       //-
    SysCtrlRegs.PCLKCR3.bit.CPUTIMER1ENCLK=0;   //-
    SysCtrlRegs.PCLKCR3.bit.CPUTIMER2ENCLK=0;   //-
    SysCtrlRegs.PCLKCR0.bit.TBCLKSYNC = 1;      // ENABLE PERIHPERAL CLOCKS AND SYNC THEM
    EPwm1Regs.TBCTL.bit.SWFSYNC = 1;            // FIRE PWM SYNCHRONIZATION PULSE

    DELAY_US(13);

    BoostTimeKeeper.ETSEL.bit.SOCAEN = 1;   
    BoostTimeKeeper.CMPA.half.CMPA = 185; // SYNC...
    BoostTimeKeeper.CMPB = 185;
    InverterTimekeeper.CMPA.half.CMPA = 1250;

    EDIS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void disableTrans(void)
{
    PwmInverter_1.CMPA.half.CMPA = INVERTER_PRD + 1;
    PwmInverter_1.CMPB = 0;

    PwmInverter_2.CMPA.half.CMPA = INVERTER_PRD + 1;
    PwmInverter_2.CMPB = 0;

    PwmInverter_3.CMPA.half.CMPA = INVERTER_PRD + 1;
    PwmInverter_3.CMPB = 0;

    PwmBoost.CMPA.half.CMPA = BOOST_PRD;
    PwmBoost.CMPB = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DSP_INIT_CONFIG()
{
    memcpy(&RamfuncsRunStart, &RamfuncsLoadStart, &RamfuncsLoadEnd - &RamfuncsLoadStart);
    InitFlash();                        // INIT FLASH IN RAM (SARAM - L0)
    DINT;                               //DISABLE ALL CPU INTS
    EALLOW;
        InitSysCtrl();
        SysCtrlRegs.HISPCP.bit.HSPCLK = 1;  //HSPCLK PRESCALER
        InitPieCtrl();                      //INIT PIE
        EnableInterrupts();
        IER = 0x0000;                       //DISABLE ALL CPU INTS
        IFR = 0x0000;                       //CLEAR ALL CPU INT FLAGS
    EDIS;
    
    AdcPowerUpAndConfig();
    InitPieVectTable();
    ConfigDma(&AdcMirror.ADCRESULT0, &DMABuff[0]); // Custom function for DMA configuration
    ConfigGpio();
    
    EALLOW;
        MapPieVectTable();
        SysCtrlRegs.PCLKCR0.bit.TBCLKSYNC = 0;
    EDIS;

    // Config(      Handle       ,Cntmod,Maxcnt       ,phaEn,       pha,Sysdiv,Hsdiv,Intsel,Inen,Indelay,phsdir)
    ConfigPwm(&BoostTimeKeeper, 2, BOOST_PRD, 1, TKBphs,  0, 0, 4, 1, 2, 1);
    ConfigPwm(&PwmBoost, 2, BOOST_PRD, 1, PBphs,  0, 0, 0, 0, 0, 1);
    ConfigPwm(&InverterTimekeeper, 2, INVERTER_TIMEKEEPER_PRD, 1, TKIphs,  0, 0, 4, 1, 1,1); //phsdir = 0 za 80khz , 1 za 40khz
    ConfigPwm(&PwmInverter_1, 2, INVERTER_PRD, 1, PIphs,  0, 0, 1, 0, 0, 0);
    ConfigPwm(&PwmInverter_2, 2, INVERTER_PRD, 1, PIphs,  0, 0, 0, 0, 0, 0);
    ConfigPwm(&PwmInverter_3, 2, INVERTER_PRD, 1, PIphs,  0, 0, 0, 0, 0, 0);

    ConfigPwmActions(&BoostTimeKeeper, 2, 2, 2, 2, 2, 2, 2, 2);
    ConfigPwmActions(&InverterTimekeeper, 2, 2, 2, 2, 2, 2, 2, 2);
    ConfigPwmActions(&PwmInverter_1, 0, 1, 0, 2, 0, 2, 0, 1);
    ConfigPwmActions(&PwmInverter_2, 0, 1, 0, 2, 0, 2, 0, 1);
    ConfigPwmActions(&PwmInverter_3, 0, 1, 0, 2, 0, 2, 0, 1);

    ConfigTripZone(&PwmInverter_1, 1);
    ConfigTripZone(&PwmInverter_2, 0);
    ConfigTripZone(&PwmInverter_3, 0);
    ConfigTripZone(&InverterTimekeeper, 0);

    SysCtrlRegs.LPMCR0.bit.QUALSTDBY = 10; // QUALIFIER FOR LP MODE

    IER |= M_INT1;  // Enable processor reaction for Group1 interrupts - ADC
    IER |= M_INT2;  // Enable processor reaction for Group2 interrupts - TZ
    IER |= M_INT3;  // Enable processor reaction for Group3 interrupts - EPWM
    IER |= M_INT7;  // Enable processor reaction for Group7 interrupts - DMA
    IER |= M_INT9;  // Enable processor reaction for Group9 interrupts - SCI
    IER |= M_INT13; //Enable processor reaction for Group13 interrupts - CpuTimer1

    EINT;           // ENABLE ALL CPU INTERRUPTS
    ERTM;           // ENABLE REAL TIME DEBBUGING INTERRUPTS
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void forceOneShotPwm()
{
    PwmInverter_1.AQSFRC.bit.OTSFA = 1;
    PwmInverter_1.AQSFRC.bit.OTSFB = 1;

    PwmInverter_2.AQSFRC.bit.OTSFA = 1;
    PwmInverter_2.AQSFRC.bit.OTSFB = 1;

    PwmInverter_3.AQSFRC.bit.OTSFA = 1;
    PwmInverter_3.AQSFRC.bit.OTSFB = 1;

    PwmBoost.AQSFRC.bit.OTSFA = 1;
    PwmBoost.AQSFRC.bit.OTSFB = 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConfigBPF(void)
{

    //POSTUPAK SINTEZE BANDPASS FILTRA JE IZ RADA:
/*
        Vol. 2 No. 3 Sep. 1998 JOURNAL OF SHANGHAI UNIVERSITY
        A Cascade Approach to the Design of IIR Bandpass Filter
        Wan Wanggen Yu Xiaoqing Yuan Jingxian Li Yonghua
        (School of Communication and Information Engineering)
*/

    delta_w_kask=2*PI*DELTA_W_BPF;              //OVA DVA PARAMETRA DEFINISU SIRINU PROPUSNOG OPSEGA BANDPASS FILTRA
    delta_w_stop_kask=2*PI*2*DELTA_W_BPF;         //OVA DVA PARAMETRA DEFINISU SIRINU PROPUSNOG OPSEGA BANDPASS FILTRA

    ks=delta_w_stop_kask/delta_w_kask;
    N_kask=4;
    qs=2*PI*F0_FILTER/delta_w_kask;
    q=qs*0.435;                         //ovde ima neko stepenovanje pa sam sam izracunao
    delta_w_pojedinacnog=2*PI*F0_FILTER/q;
    delta_w_rel=delta_w_pojedinacnog/(2*PI*fs);
    w0_rel=2*PI*F0_FILTER/fs;
    alfa=delta_w_rel*(2+delta_w_rel);
    beta=2*cos(w0_rel);
    one_minus_alfa_beta=(1-alfa)*beta;
    two_alfa_minus=2*alfa-1;

    //pocetna inicijalizacija
    y1_bef=0;
    y2_bef=0;
    y3_bef=0;
    y4_bef=0;

    y1_befbef=0;
    y2_befbef=0;
    y3_befbef=0;
    y4_befbef=0;

    sum_bef1=0;
    sum_bef2=0;
    sum_bef3=0;
    sum_bef4=0;

    sum_befbef1=0;
    sum_befbef2=0;
    sum_befbef3=0;
    sum_befbef4=0;

    kk=0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConfigPLL(void)
{
    Ts_pll= (float)(1/fs);// (float)6/(SAMPLING_FREQ);
    theta_pll=0;
    e_filt=0;
    e_filt_preth=0;
    e_filt_bef=0;
    e_filt_befbef=0;
    e_bef=0;
    e_befbef=0;
    e_bef_demod=0;
    e_befbef_demod=0;
    e_filt_demod=0;
    e_filt_bef_demod=0;
    e_filt_befbef_demod=0;
    cos_table_index=0;
    j=0;

    error_obs=0;
    up_obs=0;
    ui_obs=0;
    ui_obs_bpf=0;

    f_obs=F0_FILTER;
    f_obs_bpf=F0_FILTER;
    kp_obs=KP_OBS;
    ki_obs=KI_OBS;
    ki_obs_bpf=KI_OBS_BPF;
    ki_obs_zajedno=ki_obs*Ts_pll; //integralno dejstvo radi samo sa trenutnom greskom
    ki_obs_bpf_zajedno=ki_obs_bpf*Ts_pll;

    w_0_notch_rel_const=2*PI*2/fs;
    w_0_notch_rel=w_0_notch_rel_const*F0_FILTER;
    b1_notch=-2*cos(w_0_notch_rel);
    r=R_NOTCH;
    a1_notch=-2*r*cos(w_0_notch_rel);
    a2_notch=r*r;

	//KOEFICIJENTI KOJI DAJU POJACANJA PETLJE KADA SE PODELE SA POLOVINOM AMPLITUDE ULAZNOG SIGNALA U PLL
    ki_pll_coeff=WN_PLL*WN_PLL/(2*PI);
    kp_pll_coeff=2*KSI_PLL*WN_PLL/(2*PI);

    ki_pll=ki_pll_coeff/SCALE_FACTOR;
    kp_pll=kp_pll_coeff/SCALE_FACTOR;
    ki_pll_zajedno=ki_pll*0.5*Ts_pll;

    f_pll=F0_FILTER;
    ui_pll=F0_FILTER;
    two_pi=2*PI;

    theta_increment=2*PI/(fs);

    error_demod_obs=0;
    RSH_demod_obs=0;
    ki_demod_zajedno=1/fs;
}
