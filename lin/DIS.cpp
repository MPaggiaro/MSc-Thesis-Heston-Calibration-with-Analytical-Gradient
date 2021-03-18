#include <levmar/levmar.h>
#include <vdt/vdtMath.h>
#include "Faddeeva.hh"
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Pade approximant
#define FAST_EXP vdt::fast_exp
#define FAST_LOG vdt::fast_log

// Non-approximant
//#define FAST_EXP exp
//#define FAST_LOG log

typedef std::vector<double> VD;
//static VD x64 = {0.0243502926634244325089558, 0.0729931217877990394495429,
//                 0.1214628192961205544703765, 0.1696444204239928180373136,
//                 0.2174236437400070841496487, 0.2646871622087674163739642,
//                 0.3113228719902109561575127, 0.3572201583376681159504426,
//                 0.4022701579639916036957668, 0.4463660172534640879849477,
//                 0.4894031457070529574785263, 0.5312794640198945456580139,
//                 0.5718956462026340342838781, 0.6111553551723932502488530,
//                 0.6489654712546573398577612, 0.6852363130542332425635584,
//                 0.7198818501716108268489402, 0.7528199072605318966118638,
//                 0.7839723589433414076102205, 0.8132653151227975597419233,
//                 0.8406292962525803627516915, 0.8659993981540928197607834,
//                 0.8893154459951141058534040, 0.9105221370785028057563807,
//                 0.9295691721319395758214902, 0.9464113748584028160624815,
//                 0.9610087996520537189186141, 0.9733268277899109637418535,
//                 0.9833362538846259569312993, 0.9910133714767443207393824,
//                 0.9963401167719552793469245, 0.9993050417357721394569056};
//static VD w64 = {0.0486909570091397203833654, 0.0485754674415034269347991,
//                 0.0483447622348029571697695, 0.0479993885964583077281262,
//                 0.0475401657148303086622822, 0.0469681828162100173253263,
//                 0.0462847965813144172959532, 0.0454916279274181444797710,
//                 0.0445905581637565630601347, 0.0435837245293234533768279,
//                 0.0424735151236535890073398, 0.0412625632426235286101563,
//                 0.0399537411327203413866569, 0.0385501531786156291289625,
//                 0.0370551285402400460404151, 0.0354722132568823838106931,
//                 0.0338051618371416093915655, 0.0320579283548515535854675,
//                 0.0302346570724024788679741, 0.0283396726142594832275113,
//                 0.0263774697150546586716918, 0.0243527025687108733381776,
//                 0.0222701738083832541592983, 0.0201348231535302093723403,
//                 0.0179517157756973430850453, 0.0157260304760247193219660,
//                 0.0134630478967186425980608, 0.0111681394601311288185905,
//                 0.0088467598263639477230309, 0.0065044579689783628561174,
//                 0.0041470332605624676352875, 0.0017832807216964329472961};
//
static VD x128 = {0.0122236989606157641980521,0.0366637909687334933302153,0.0610819696041395681037870,0.0854636405045154986364980,0.1097942311276437466729747,0.1340591994611877851175753,0.1582440427142249339974755,0.1823343059853371824103826,0.2063155909020792171540580,0.2301735642266599864109866,0.2538939664226943208556180,0.2774626201779044028062316,0.3008654388776772026671541,0.3240884350244133751832523,0.3471177285976355084261628,0.3699395553498590266165917,0.3925402750332674427356482,0.4149063795522750154922739,0.4370245010371041629370429,0.4588814198335521954490891,0.4804640724041720258582757,0.5017595591361444642896063,0.5227551520511754784539479,0.5434383024128103634441936,0.5637966482266180839144308,0.5838180216287630895500389,0.6034904561585486242035732,0.6228021939105849107615396,0.6417416925623075571535249,0.6602976322726460521059468,0.6784589224477192593677557,0.6962147083695143323850866,0.7135543776835874133438599,0.7304675667419088064717369,0.7469441667970619811698824,0.7629743300440947227797691,0.7785484755064119668504941,0.7936572947621932902433329,0.8082917575079136601196422,0.8224431169556438424645942,0.8361029150609068471168753,0.8492629875779689691636001,0.8619154689395484605906323,0.8740527969580317986954180,0.8856677173453972174082924,0.8967532880491581843864474,0.9073028834017568139214859,0.9173101980809605370364836,0.9267692508789478433346245,0.9356743882779163757831268,0.9440202878302201821211114,0.9518019613412643862177963,0.9590147578536999280989185,0.9656543664319652686458290,0.9717168187471365809043384,0.9771984914639073871653744,0.9820961084357185360247656,0.9864067427245862088712355,0.9901278184917343833379303,0.9932571129002129353034372,0.9957927585349811868641612,0.9977332486255140198821574,0.9990774599773758950119878,0.9998248879471319144736081};
static VD w128 = {0.0244461801962625182113259,0.0244315690978500450548486,0.0244023556338495820932980,0.0243585572646906258532685,0.0243002001679718653234426,0.0242273192228152481200933,0.0241399579890192849977167,0.0240381686810240526375873,0.0239220121367034556724504,0.0237915577810034006387807,0.0236468835844476151436514,0.0234880760165359131530253,0.0233152299940627601224157,0.0231284488243870278792979,0.0229278441436868469204110,0.0227135358502364613097126,0.0224856520327449668718246,0.0222443288937997651046291,0.0219897106684604914341221,0.0217219495380520753752610,0.0214412055392084601371119,0.0211476464682213485370195,0.0208414477807511491135839,0.0205227924869600694322850,0.0201918710421300411806732,0.0198488812328308622199444,0.0194940280587066028230219,0.0191275236099509454865185,0.0187495869405447086509195,0.0183604439373313432212893,0.0179603271850086859401969,0.0175494758271177046487069,0.0171281354231113768306810,0.0166965578015892045890915,0.0162550009097851870516575,0.0158037286593993468589656,0.0153430107688651440859909,0.0148731226021473142523855,0.0143943450041668461768239,0.0139069641329519852442880,0.0134112712886163323144890,0.0129075627392673472204428,0.0123961395439509229688217,0.0118773073727402795758911,0.0113513763240804166932817,0.0108186607395030762476596,0.0102794790158321571332153,0.0097341534150068058635483,0.0091830098716608743344787,0.0086263777986167497049788,0.0080645898904860579729286,0.0074979819256347286876720,0.0069268925668988135634267,0.0063516631617071887872143,0.0057726375428656985893346,0.0051901618326763302050708,0.0046045842567029551182905,0.0040162549837386423131943,0.0034255260409102157743378,0.0028327514714579910952857,0.0022382884309626187436221,0.0016425030186690295387909,0.0010458126793403487793129,0.0004493809602920903763943};

typedef std::complex<double> CD;
CD one(1.0, 0.0), zero(0.0, 0.0), two(2.0, 0.0), i(0.0, 1.0);
const double pi = 4.0 * atan(1.0), lb = 0.0, ub = 200, mid = 0.5 * (ub + lb),
             halfRange = 0.5 * (ub - lb);

struct GLsetting {
    int nGrid;
    VD *abs;
    VD *weights;
};

//static GLsetting gl = {64, &x64, &w64};
static GLsetting gl = {128, &x128, &w128};

struct CFPriceData {
    CD CFPrice;
    CD iu;
    CD interU;
    CD xi;
    CD d;
    CD interD;
    CD A1;
    CD A2;
    CD A;
    CD D;
};

struct CFvolData {
    CD CFvol;
    CD G;
    CD F;
};

struct intVIXdata {
    double VIXint;
    CD rawVIXint;
    double atauBar;
    CD G;
    CD F;
    CD inputU;
    CD iu;
    int start;
    int end;
};

struct intI {
    double start;
    double end;
    double inte;
};

struct mktPara {
    double S;
    double r;
    VD SPX_T;
    VD SPX_K;
    VD VIX_T;
    VD VIX_K;
    VD mktPrices;
    double tbar;
    std::vector<intI> toCali;
};

struct modelPara {
    double k;
    double vbar;
    double v0;
    double rho;
    double sigma;
    std::vector<intI> disArr;
};

CFPriceData CFprice(CD u, modelPara p, double tau, double S, double r);
double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r);
VD SPXprice(modelPara p, VD tau, double S, VD K, double r, int n);
VD gradSPXintgrand(double u, modelPara p, double tau, double K, double S,
                   double r);
VD gradientSPXprice(modelPara p, double S, double r, int n, VD tau, VD K);

CFvolData CFvol(CD u, modelPara p, double tau);
intVIXdata VIXintegrand(CD u, modelPara p, double tau, double K, double tbar);
VD VIXprice(modelPara p, VD tau, double tbar, VD K, double r, int n);
VD gradVIXintegrand(CD u, modelPara p, double tau, double K, double tbar);
VD gradientVIXprice(modelPara p, double r, int n, VD tau, VD K, double tbar);

void objFunc(double *p, double *x, int m, int n, void *data);
void JacFunc(double *p, double *jac, int m, int n, void *data);
// For comparing
void showSPXcallPrices(modelPara mp, VD tarr, double S, VD karr, double r,
                       int n);
void printCFvol(modelPara p, double tau, int n);
void printIntegrandVIXoption(modelPara p, double tau, double K, double tbar,
                             int n);
void printVIXcalls(modelPara p, VD tau, double tbar, VD K, double r, int n);
void printSPXgradient(modelPara p, double S, double r, int n, VD tau, VD K);
void printVIXgradient(modelPara p, double tbar, double r, int n, VD tau, VD K);

int main() {
    // These are manully calculated according to the SPXtarr and VIXtarr
    VD allT = {0.0,        0.01917808, 0.03835616, 0.05753425, 0.08219178,
               0.10136986, 0.12054795, 0.13972603, 0.16438356, 0.24657534,
               0.32876712, 0.4109589,  0.49315068, 0.73972603, 0.98630137};

    // These are uniform random (0, 1.4) generated
    VD disStar = {1.46371706e-04, 4.31183274e-04, 7.89958945e-04, 2.91974407e-04,
       8.57045896e-05, 5.48371267e-04, 4.84620007e-04, 3.51856480e-04,
       7.00281831e-04, 7.14263885e-04, 7.51713276e-04, 6.58928904e-04,
       4.95814626e-04, 6.67154021e-04};

    std::vector<intI> disArr;

    for (int j = 0; j < 14; j++) {
        disArr.push_back({allT[j], allT[j + 1], disStar[j]});
    }

    modelPara mp = {3.0, 0.1, 0.08, -0.8, 0.25, disArr};
    double S = 1.0;
    double r = 0.02;
    double tbar = 30 / 365.0;

    // Maturities for VIX options: 7, 14, 21, 30, 60, 90, (divided by 365)
    // BSdelta: 0.9, 0.75, 0.5, 0.25, 0.1
    const VD VIXkarr = {
        0.07660333466800909, 0.07823036262505452, 0.08007866879151475,
        0.08197064388101638, 0.08371167160297435, 0.07528231789011931,
        0.07755350918528185, 0.08015741494276396, 0.08284874840358378,
        0.0853482112450969,  0.07430112345712658, 0.077055734551228,
        0.08023623852982013, 0.08354801898792094, 0.08664544590808211,
        0.0732872469748735,  0.07654634796281677, 0.08033769704753065,
        0.08431683207193616, 0.08806642128491488, 0.07084911309063555,
        0.07534553012964565, 0.08067681958626018, 0.08638533974549638,
        0.0918677586016848,  0.06910110193769552, 0.07450988995253548,
        0.08101737363349071, 0.08809320258894406, 0.0949885695945464};

    const VD VIXtarr = {
        0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.03835616,
        0.03835616, 0.03835616, 0.03835616, 0.03835616, 0.05753425, 0.05753425,
        0.05753425, 0.05753425, 0.05753425, 0.08219178, 0.08219178, 0.08219178,
        0.08219178, 0.08219178, 0.16438356, 0.16438356, 0.16438356, 0.16438356,
        0.16438356, 0.24657534, 0.24657534, 0.24657534, 0.24657534, 0.24657534};

    const VD SPXkarr = {0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699,
                        0.6137, 0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358,
                        0.9175, 0.9025, 1.0427, 1.0463, 1.0499, 1.0530, 1.0562,
                        1.0593, 1.0663, 1.0766, 1.2287, 1.2399, 1.2485, 1.2659,
                        1.2646, 1.2715, 1.2859, 1.3046, 1.3939, 1.4102, 1.4291,
                        1.4456, 1.4603, 1.4736, 1.5005, 1.5328};

    const VD SPXtarr = {
        0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
        0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
        0.4109589,  0.49315068, 0.73972603, 0.98630137, 0.08219178, 0.16438356,
        0.24657534, 0.32876712, 0.4109589,  0.49315068, 0.73972603, 0.98630137,
        0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
        0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
        0.4109589,  0.49315068, 0.73972603, 0.98630137};

    //
    VD SPXprices = SPXprice(mp, SPXtarr, S, SPXkarr, r, (int)SPXkarr.size());
    VD VIXprices = VIXprice(mp, VIXtarr, tbar, VIXkarr, r, (int)VIXkarr.size());
    VD optPrices;
    optPrices.reserve(SPXprices.size() + VIXprices.size());
    optPrices.insert(optPrices.end(), SPXprices.begin(), SPXprices.end());
    optPrices.insert(optPrices.end(), VIXprices.begin(), VIXprices.end());

    std::vector<intI> toCali;

    for (int j = 0; j < 14; j++) {
        toCali.push_back({allT[j], allT[j + 1], 0.0});
    }
    mktPara marP = {S, r, SPXtarr, SPXkarr, VIXtarr, VIXkarr, optPrices,
    tbar, toCali};
    //mktPara marP = {S, r, SPXtarr, SPXkarr, VIXtarr, VIXkarr, optPrices,
    //tbar, disArr};


    VD disInitial = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    double p[19];
    // Initial 1
    p[0] = 1.2000;
    p[1] = 0.20000;
    p[2] = 0.2000;
    p[3] = -0.6000;
    p[4] = 0.3000;

    // Initial 2
    //p[0] = 1.2000;
    //p[1] = 0.30000;
    //p[2] = 0.2500;
    //p[3] = -0.4000;
    //p[4] = 0.2000;

    // Initial 3
    //p[0] = 1.4000;
    //p[1] = 0.40000;
    //p[2] = 0.4000;
    //p[3] = -0.7000;
    //p[4] = 0.2000;
    for(int fill = 0; fill < 14; fill++){
        p[fill+5] = disInitial[fill] * 0.0001;
    }



    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0] = LM_INIT_MU * 100;
    // stopping thresholds for
    opts[1] = 1E-10;          // ||J^T e||_inf
    opts[2] = 1E-300;          // ||Dp||_2
    opts[3] = 1E-10;          // ||e||_2
    opts[4] = LM_DIFF_DELTA;  // finite difference if used

    std::cout <<
    "\r-------- -------- -------- Heston Model Calibrator-------- -------- --------"
    <<std::endl;
    
    std::cout <<
    "Parameters:" << "\t      kappa"<<"\t     vinf"<< "\t       v0"<<"\t      rho"
    <<"\t       vov"<<"\t       I1"<< "\t       I2"<< "\t       I3"<< "\t       I4"
    <<"\t       I5"<<"\t       I6"<<"\t       I7"<<"\t       I8"<<"\t       I9"
    <<"\t       I10"<<"\t       I11"<<"\t       I12"<<"\t       I13"<<"\t       I14"
    <<std::endl;

    std::cout << "\r Initial point:" << "\t"  << std::scientific
    << std::setprecision(8) << p[0]<< "\t" << p[1]<< "\t"<< p[2]<< "\t"<< p[3]
    << "\t"<< p[4] << "\t" << p[5]<< "\t" << p[6]<< "\t"<< p[7]<< "\t"<< p[8]
    << "\t"<< p[9] << "\t" << p[10]<< "\t" << p[11]<< "\t"<< p[12]<< "\t"<< p[13]
    << "\t"<< p[14] << "\t" << p[15]<< "\t" << p[16]<< "\t"<< p[17]<< "\t"<< p[18]
    << std::endl;

    double err[optPrices.size()];
    dlevmar_chkjac(objFunc, JacFunc, p, 5 + 14, optPrices.size(), (void *)&marP, err);
    for (int cer = 0; cer < 30; cer++){
        std::cerr << err[cer] << ",";
    }
    std::cerr << std::endl;

    double start_s = clock();
    dlevmar_der(objFunc, JacFunc, p, NULL, 19, (int)optPrices.size(), 30000,
    opts, info, NULL, NULL, (void *) &marP);
    //dlevmar_dif(objFunc, p, NULL, 19, (int)optPrices.size(), 30000,
    //opts, info, NULL, NULL, (void *) &marP);
    double stop_s = clock();

    std::cout << "Optimum found:" << std::scientific << std::setprecision(8)
    << "\t"<< p[0]<< "\t" << p[1]<< "\t"<< p[2]<< "\t"<< p[3]<< "\t"<< p[4]
    << "\t" << p[5]<< "\t" << p[6]<< "\t"<< p[7]<< "\t"<< p[8]
    << "\t"<< p[9] << "\t" << p[10]<< "\t" << p[11]<< "\t"<< p[12]<< "\t"<< p[13]
    << "\t"<< p[14] << "\t" << p[15]<< "\t" << p[16]<< "\t"<< p[17]<< "\t"<< p[18]
    << std::endl;

    std::cout << "Real optimum:" << "\t" << mp.k<<"\t"<< mp.vbar
    << "\t"<< mp.v0<< "\t"<< mp.rho<< "\t"<< mp.sigma << disStar[0]
    << "\t"<< disStar[1]<< "\t"<< disStar[2]<< "\t"<< disStar[3]
    << "\t"<< disStar[4]<< "\t"<< disStar[5]<< "\t"<< disStar[6]
    << "\t"<< disStar[7]<< "\t"<< disStar[8]<< "\t"<< disStar[9]
    << "\t"<< disStar[10]<< "\t"<< disStar[11]<< "\t"<< disStar[12]
    << "\t"<< disStar[13] << std::endl;


    if (int(info[6]) == 6) {
        std::cout << "\r Solved: stopped by small ||e||_2 = " << info[1] << " < "
             << opts[3] << std::endl;
    } else if (int(info[6]) == 1) {
        std::cout << "\r Solved: stopped by small gradient J^T e = " << info[2]
             << " < " << opts[1] << std::endl;
    } else if (int(info[6]) == 2) {
        std::cout << "\r Solved: stopped by small change Dp = " << info[3] << " < "
             << opts[2] << std::endl;
    } else if (int(info[6]) == 3) {
        std::cout << "\r Unsolved: stopped by itmax " << std::endl;
    } else if (int(info[6]) == 4) {
        std::cout << "\r Unsolved: singular matrix. Restart from current p with "
                "increased mu"
             << std::endl;
    } else if (int(info[6]) == 5) {
        std::cout << "\r Unsolved: no further error reduction is possible. Restart "
                "with increased mu"
             << std::endl;
    } else if (int(info[6]) == 7) {
        std::cout << "\r Unsolved: stopped by invalid values, user error" << std::endl;
    }

    std::cout << "\r-------- -------- -------- Computational cost -------- -------- "
            "--------"
         << std::endl;
    std::cout << "\r          Time cost: "
         << double(stop_s - start_s) / CLOCKS_PER_SEC << " seconds " << std::endl;
    std::cout << "         Iterations: " << int(info[5]) << std::endl;
    std::cout << "         pv  Evalue: " << int(info[7]) << std::endl;
    std::cout << "         Jac Evalue: " << int(info[8]) << std::endl;
    std::cout << "# of lin sys solved: " << int(info[9])
         << std::endl;  // The attempts to reduce error
    std::cout << "\r-------- -------- -------- Residuals -------- -------- --------"
         << std::endl;
    std::cout << " \r            ||e0||_2: " << info[0] << std::endl;
    std::cout << "           ||e*||_2: " << info[1] << std::endl;
    std::cout << "          ||J'e||_inf: " << info[2] << std::endl;
    std::cout << "           ||Dp||_2: " << info[3] << std::endl;

    //showSPXcallPrices(mp, SPXtarr, S, SPXkarr, r, (int)SPXkarr.size());
    //printSPXgradient(mp, S, r, (int)SPXkarr.size(), SPXtarr, SPXkarr);
    //printCFvol(mp, VIXtarr[1], gl.nGrid>>1);
    //printIntegrandVIXoption(mp, VIXtarr[5], VIXkarr[5], tbar, gl.nGrid>>1);
    printVIXcalls(mp, VIXtarr, tbar, VIXkarr, r, (int)VIXkarr.size());
    //printVIXgradient(mp, tbar, r, (int)VIXkarr.size(), VIXtarr, VIXkarr);
    //CD input = mid + x64[0] * halfRange + i;
    //gradVIXintegrand(input, mp, VIXtarr[1], VIXkarr[1], tbar);
    //VIXintegrand(input, mp, VIXtarr[25], VIXkarr[25], tbar);
    //
    //
    //for (int kk = 0; kk < 30; kk++){
    //    VIXintegrand(x128[10] + i, mp, VIXtarr[kk], VIXkarr[kk], tbar);
    //}
}

// Functions to keep
void objFunc(double *p, double *x, int m, int n, void *data) {
    int k;
    (void)m;
    mktPara *mktp;
    mktp = (struct mktPara *)data;
    std::vector<intI> toCali = mktp->toCali;

    for(int count = 0; count < 14; count++){
        toCali[count].inte = p[count+5];
    }
    //for(int print = 0; print < 14; print++){
    //    std::cout << "In objFunc " << print << " " << toCali[print].inte << std::endl;
    //}
    modelPara molp = {p[0], p[1], p[2], p[3], p[4], toCali};
    VD SPXprices = SPXprice(molp, mktp->SPX_T, mktp->S, mktp->SPX_K, mktp->r,
                            (int)mktp->SPX_T.size());
    VD VIXprices = VIXprice(molp, mktp->VIX_T, mktp->S, mktp->VIX_K, mktp->r,
                            (int)mktp->VIX_K.size());
    double Nspx2root = sqrt(2 * SPXprices.size());
    double Nvix2root = sqrt(2 * VIXprices.size());
    for (k = 0; k < (int)SPXprices.size(); k++) {
        x[k] = (SPXprices[k] - mktp->mktPrices[k]) /
               (mktp->mktPrices[k] * Nspx2root);
    }
    for (int j = k; j < n; j++) {
        x[j] = (VIXprices[j - k] - mktp->mktPrices[j]) /
               (mktp->mktPrices[j] * Nvix2root);
    }
}

void JacFunc(double *p, double *jac, int m, int n, void *data) {
    mktPara *mktp;
    mktp = (struct mktPara *)data;
    std::vector<intI> toCali = mktp->toCali;

    for(int count = 0; count < 14; count++){
        toCali[count].inte = p[count+5];
    }
    //for(int print = 0; print < 14; print++){
    //    std::cout << "In JacFunc " << print << " " << toCali[print].inte << std::endl;
    //}
    modelPara molp = {p[0], p[1], p[2], p[3], p[4], toCali};

    VD SPXjac =
        gradientSPXprice(molp, mktp->S, mktp->r, (int)mktp->SPX_T.size(),
                         mktp->SPX_T, mktp->SPX_K);
    VD VIXjac = gradientVIXprice(molp, mktp->r, (int)mktp->VIX_T.size(),
                                 mktp->VIX_T, mktp->VIX_K, mktp->tbar);
    double Nspx2root = sqrt(2 * mktp->SPX_T.size());
    double Nvix2root = sqrt(2 * mktp->VIX_T.size());
    int k;
    for (k = 0; k < (int)mktp->SPX_T.size(); k++) {
        for (int j = 0; j < m; j++) {
            jac[k * m + j] =
                SPXjac[k * m + j] / (mktp->mktPrices[k] * Nspx2root);
        }
    }
    for (int l = k; l < n; l++) {
        for (int j = 0; j < m; j++) {
            jac[l * m + j] =
                VIXjac[(l - k) * m + j] / (mktp->mktPrices[l] * Nvix2root);
        }
    }
}

// CF of price, stock option pircing, and gradient of stock option price
CFPriceData CFprice(CD u, modelPara p, double tau, double S, double r) {
    double var = pow(p.sigma, 2);

    CD iu = i * u;
    CD interU = pow(u, 2) + iu;
    CD xi = p.k - p.sigma * p.rho * iu;
    CD d = sqrt(pow(xi, 2) + var * interU);
    CD interD = d * tau * 0.5;
    CD A1 = interU * sinh(interD);
    CD A2 = (d * cosh(interD) + xi * sinh(interD)) / p.v0;
    CD A = A1 / A2;
    CD D = log(d / (p.v0 * A2)) + p.k * tau * 0.5;
    double tmp = p.k * p.vbar / p.sigma;

    CD CF = exp(iu * (FAST_LOG(S) + r * tau) - tmp * p.rho * tau * iu - A +
                2 * tmp * D / p.sigma);

    CFPriceData ret = {CF, iu, interU, xi, d, interD, A1, A2, A, D};
    return ret;
}

double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r) {
    CD iu = i * u;
    CD inputU = u - i * 0.5;
    double x = FAST_LOG(S);
    double rT = r * tau;
    double kappa = x - FAST_LOG(K) + rT;
    double Iphi = 0.0;
    int count = 0;
    do {
        Iphi += p.disArr.at(count).inte;
        count++;
    } while (count < 14 && p.disArr.at(count).end <= tau);

    CD integrand1 = exp(iu * kappa - i * inputU * (x + rT));
    CFPriceData tmp = CFprice(inputU, p, tau, S, r);
    CD CFpri = tmp.CFPrice;
    double tmpU = (pow(u, 2) + 0.25);
    double integrand2 = exp(-tmpU * Iphi);
    double SPXint = real(integrand1 * CFpri) * integrand2 / tmpU;

    return SPXint;
}

VD SPXprice(modelPara p, VD tau, double S, VD K, double r,
            int n) {  // tau and K can be passed by reference, but we will see,
                      // may be I'll make them const or static
    int nGrid = gl.nGrid;

    double up_u, down_u, upInt, downInt, strike, T, rT, glCollect, glInt,
        SPXcall;
    nGrid = nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD SPXs;
    SPXs.reserve(n);

    for (int j = 0; j < n; j++) {
        strike = K[j];
        T = tau[j];
        rT = r * T;
        glCollect = 0.0;

        for (int count = 0; count < nGrid; count++) {
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = SPXintegrand(up_u, p, T, strike, S, r);
            downInt = SPXintegrand(down_u, p, T, strike, S, r);
            glCollect += w[count] * (upInt + downInt);
        }

        glInt = halfRange * glCollect;
        SPXcall = S - sqrt(strike * S) * FAST_EXP(-rT * 0.5) / pi * glInt;
        SPXs.push_back(SPXcall);
    }

    return SPXs;  // Here can return an adress of the VD.
}

VD gradSPXintgrand(double u, modelPara p, double tau, double K, double S,
                   double r) {
    VD gradSPXint;
    gradSPXint.reserve(5);
    CD inputU = u - i * 0.5;
    double var = pow(p.sigma, 2);
    CFPriceData CFP = CFprice(inputU, p, tau, S, r);
    CD B = exp(CFP.D);
    CD tmpBase = p.k * tau * CFP.iu / p.sigma;
    double tmp1 = 2 * p.k * p.vbar / var;
    CD dOverA2 = CFP.d / CFP.A2;
    CD interXi = 2.0 + tau * CFP.xi;
    CD revA2 = 1.0 / CFP.A2;
    CD AoverA2 = CFP.A / CFP.A2;
    CD coshInterD = cosh(CFP.interD);

    // Partials for rho
    CD d_rho = -CFP.xi * p.sigma * CFP.iu / CFP.d;
    CD A1_rho = -CFP.iu * CFP.interU * tau * CFP.xi * p.sigma / (2.0 * CFP.d) *
                coshInterD;
    CD A2_rho = -p.sigma * CFP.iu * interXi / (2.0 * CFP.d * p.v0) *
                (CFP.xi * coshInterD + CFP.d * sinh(CFP.interD));
    CD A_rho = revA2 * A1_rho - AoverA2 * A2_rho;
    CD B_rho = exp(p.k * tau * 0.5) / p.v0 *
               (revA2 * d_rho - dOverA2 / CFP.A2 * A2_rho);

    // Partials for k
    CD B_k = i / (p.sigma * inputU) * B_rho + B * tau * 0.5;

    // Partials for sigma
    CD d_sigma = (p.rho / p.sigma - 1.0 / CFP.xi) * d_rho +
                 p.sigma * pow(inputU, 2) / CFP.d;
    CD A1_sigma = CFP.interU * tau * d_sigma * coshInterD / 2.0;
    CD A2_sigma = p.rho / p.sigma * A2_rho -
                  interXi / (p.v0 * tau * CFP.xi * CFP.iu) * A1_rho +
                  p.sigma * tau * CFP.A1 / (2.0 * p.v0);
    CD A_sigma = revA2 * A1_sigma - AoverA2 * A2_sigma;

    // Gradient of CF price
    CD v0Par = (-CFP.A / p.v0) * CFP.CFPrice;
    CD vbarPar = (tmp1 / p.vbar * CFP.D - tmpBase * p.rho) * CFP.CFPrice;
    CD rhoPar = (-A_rho + tmp1 / CFP.d * (d_rho - dOverA2 * A2_rho) -
                 tmpBase * p.vbar) *
                CFP.CFPrice;
    CD kPar = (1.0 / (p.sigma * CFP.iu) * A_rho + tmp1 / p.k * CFP.D +
               tmp1 / B * B_k - tmpBase / p.k * p.vbar * p.rho) *
              CFP.CFPrice;
    CD sigmaPar = (-A_sigma - 2.0 * tmp1 / p.sigma * CFP.D +
                   tmp1 / CFP.d * (d_sigma - dOverA2 * A2_sigma) +
                   tmpBase * p.vbar * p.rho / p.sigma) *
                  CFP.CFPrice;

    double x = FAST_LOG(S);
    double rT = r * tau;
    double kappa = x - FAST_LOG(K) + rT;
    CD integrand1 = exp(i * u * kappa - i * inputU * (x + rT));
    double integrand2 = pow(u, 2) + 0.25;
    double Iphi = 0.0;
    int count = 0;
    do {
        Iphi += p.disArr.at(count).inte;
        count++;
    } while (count < 14 && p.disArr.at(count).end <= tau);
    double integrand3 = exp(-integrand2 * Iphi);

    gradSPXint.push_back(real(integrand1 * kPar) * integrand3 / integrand2);
    gradSPXint.push_back(real(integrand1 * vbarPar) * integrand3 / integrand2);
    gradSPXint.push_back(real(integrand1 * v0Par) * integrand3 / integrand2);
    gradSPXint.push_back(real(integrand1 * rhoPar) * integrand3 / integrand2);
    gradSPXint.push_back(real(integrand1 * sigmaPar) * integrand3 / integrand2);

    double IphiDer = -real(integrand1 * CFP.CFPrice) * integrand3;
    for (int new_count = 0; new_count < count; new_count++) {
        gradSPXint.push_back(IphiDer);
    }
    int to_fill = 19 - (int)gradSPXint.size();
    if (to_fill > 0) {  // this is test necessary
        for (int fill = 0; fill < to_fill; fill++) {
            gradSPXint.push_back(0.0);
        }
    }

    return gradSPXint;
}

VD gradientSPXprice(modelPara p, double S, double r, int n, VD tau, VD K) {
    int nGrid = gl.nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD gradSPX, upInt, downInt, glCollect;
    gradSPX.reserve(19 * n);
    glCollect.reserve(19);

    double strike, T, rT, up_u, down_u;
    for (int l = 0; l < n; l++) {
        strike = K[l];
        T = tau[l];
        rT = r * T;
        for (int cc = 0; cc < 19; cc++) {
            glCollect[cc] = 0.0;
        }  // There must be a better way
        for (int count = 0; count < nGrid; count++) {
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = gradSPXintgrand(up_u, p, T, strike, S, r);
            downInt = gradSPXintgrand(down_u, p, T, strike, S, r);
            for (int j = 0; j < 19; j++)
                glCollect[j] += w[count] * (upInt[j] + downInt[j]);
        }
        for (int p = 0; p < 19; p++) {
            glCollect[p] = glCollect[p] * halfRange;
            gradSPX.push_back(-sqrt(strike * S) * FAST_EXP(-rT * 0.5) / pi *
                              glCollect[p]);
        }
    }

    return gradSPX;
}

// CF of vol, VIX option pricing, and gradient of VIX option pricing
CFvolData CFvol(CD u, modelPara p, double tau) {
    double var = pow(p.sigma, 2);

    CD iu = u * i;
    double ktauH = p.k * tau * 0.5;
    CD G = cosh(ktauH) + (1.0 - var * iu / p.k) * sinh(ktauH);
    CD F = p.v0 * iu / G * FAST_EXP(-ktauH);

    CD CF = pow(FAST_EXP(ktauH) / G, 2 * p.k * p.vbar / var) * exp(F);
    CFvolData ret = {CF, G, F};

    return ret;
}

intVIXdata VIXintegrand(CD u, modelPara p, double tau, double K, double tbar) {
    CD iu = i * u;
    double k = p.k;
    double vbar = p.vbar;
    double atauBar = (1.0 - FAST_EXP(-tbar * k)) / k;
    CD inputU = -u * atauBar / tbar;
    double btauBar = vbar * (tbar - atauBar);
    CFvolData tmp = CFvol(inputU, p, tau);
    CD CFvola = tmp.CFvol;
    double Iphi = 0.0;
    int count = 0;
    while (p.disArr.at(count).start < tau && count < 13) {
        count++;
    }
    int start_count = count;
    do {
        Iphi += p.disArr.at(count).inte;
        count++;
    } while (p.disArr.at(count).end <= (tau + 30.0 / 365 + 0.00001) &&
             count < 13);
    int end_count = count - 1;

    CD part1 = exp(-iu * (btauBar + Iphi) / tbar);
    CD part2 = 1.0 - Faddeeva::erf(K * sqrt(-iu));
    CD part3 = pow(-iu, 3 / 2.0);
    CD rawInt = CFvola * part1 * part2 / part3;

    double VIXint = real(rawInt);

    intVIXdata ret = {VIXint, rawInt, atauBar,     tmp.G,    tmp.F,
                      inputU, iu,     start_count, end_count};

    return ret;
}

// For this pricing, the u is not a real but can be a complex. See smile page
// 7, equation (11). Currently, u is real which means im(u) is 0, but in the
// paper it said it should be > 0. Do I choose one im(u) and how to choose?
// So does this numerical integration still work? And what should the upper
// boundary be?
VD VIXprice(modelPara p, VD tau, double tbar, VD K, double r, int n) {
    int nGrid = gl.nGrid;

    double up_u, down_u, upInt, downInt, strike, T, discount, glCollect, glInt,
        VIXcall;
    nGrid = nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD VIXs;
    VIXs.reserve(n);

    for (int j = 0; j < n; j++) {
        strike = K[j];
        T = tau[j];
        discount = exp(-r * T);
        glCollect = 0.0;
        for (int count = 0; count < nGrid; count++) {
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = VIXintegrand(up_u + i, p, T, strike, tbar).VIXint;
            // What should be the complex part be? By adding
            // this magic number I already made the price
            // positive, but how to choose the complex part?
            downInt = VIXintegrand(down_u + i, p, T, strike, tbar).VIXint;
            glCollect += w[count] * (upInt + downInt);
        }

        glInt = halfRange * glCollect;
        VIXcall = 0.5 * discount / sqrt(pi) * glInt;
        VIXs.push_back(VIXcall);
    }
    return VIXs;
}

VD gradVIXintegrand(CD u, modelPara p, double tau, double K, double tbar) {
    intVIXdata inter = VIXintegrand(u, p, tau, K, tbar);
    double var = pow(p.sigma, 2);
    double tmp1 = 2.0 * p.k / var;
    double tmp2 = p.k * tau * 0.5;
    double tmp3 = FAST_EXP(tmp2);
    CD tmp4 = pow(inter.G, 2);
    CD iU = i * inter.inputU;
    CD iuOverTbar = inter.iu / tbar;

    // equation (3.34)
    CD G_sigma = -2.0 * p.sigma * iU / p.k * sinh(tmp2);

    CD h_v0 = inter.F / p.v0;
    CD h_vbar = tmp1 * log(tmp3 / inter.G);
    CD h_sigma = -2.0 * p.vbar / p.sigma * h_vbar -
                 tmp1 * p.vbar / inter.G * G_sigma -
                 p.v0 * iU / (tmp4 * tmp3) * G_sigma;
    CD h_k = -p.sigma / (2.0 * p.k) * h_sigma +
             p.vbar * tau * iU / (inter.G * tmp3) -
             p.v0 * inter.inputU * tau / (2.0 * p.k * tmp4) *
                 (2.0 * p.k * i + inter.inputU * var);

    // equation (3.16)
    double atauBar_k = (tbar - inter.atauBar * (p.k * tbar + 1)) / p.k;
    double btauBar_vbar = tbar - inter.atauBar;
    double btauBar_k = -p.vbar * atauBar_k;

    // equation (3.23), (3.24)
    CD G_U = (inter.G - tmp3) / inter.inputU;
    CD hPrime = h_k - u / tbar *
                          (inter.F / inter.inputU -
                           1.0 / inter.G * (tmp1 * p.vbar + inter.F) * G_U) *
                          atauBar_k;

    // equation (3.29)
    CD H_vbar = h_vbar - iuOverTbar * btauBar_vbar;
    CD H_k = hPrime - iuOverTbar * btauBar_k;
    double IphiDer = real(inter.rawVIXint * (-inter.iu) / tbar);

    VD ret;
    ret.reserve(5);

    // equation (3.28)
    ret.push_back(real(H_k * inter.rawVIXint));
    ret.push_back(real(H_vbar * inter.rawVIXint));
    ret.push_back(real(h_v0 * inter.rawVIXint));
    ret.push_back(0.0);
    ret.push_back(real(h_sigma * inter.rawVIXint));

    int count;
    for (count = 0; count < inter.start; count++) {
        ret.push_back(0.0);
    }

    for (count = inter.start; count <= inter.end; count++) {
        ret.push_back(IphiDer);
    }
    for (count = inter.end; count < 14; count++) {
        ret.push_back(0.0);
    }

    return ret;
}

VD gradientVIXprice(modelPara p, double r, int n, VD tau, VD K, double tbar) {
    int nGrid = gl.nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD gradVIX, upInt, downInt, glCollect;
    gradVIX.reserve(19 * n);
    glCollect.reserve(19);

    double strike, T, discount, up_u, down_u;
    double fixedPart = 0.5 / sqrt(pi);
    for (int l = 0; l < n; l++) {
        strike = K[l];
        T = tau[l];
        discount = exp(-r * T);
        for (int cc = 0; cc < 19; cc++) glCollect[cc] = 0.0;
        for (int count = 0; count < nGrid; count++) {
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = gradVIXintegrand(up_u + i, p, T, strike, tbar);
            downInt = gradVIXintegrand(down_u + i, p, T, strike, tbar);
            for (int j = 0; j < 19; j++)
                glCollect[j] = w[count] * (upInt[j] + downInt[j]);
        }
        for (int p = 0; p < 19; p++) {
            glCollect[p] = glCollect[p] * halfRange;
            gradVIX.push_back(fixedPart * discount * glCollect[p]);
        }
    }

    return gradVIX;
}

// Functions for comparing results
void printSPXgradient(modelPara p, double S, double r, int n, VD tau, VD K) {
    VD gradients = gradientSPXprice(p, S, r, n, tau, K);
    for (int j = 0; j < 40; j++) {
        for (int jj = 0; jj < 5; jj++){
            std::cout << std::setprecision(16);
            std::cout << "gradient" << j*5+jj << "   " << gradients[j*5+jj] << std::endl;
        }
    }
}

void printVIXgradient(modelPara p, double tbar, double r, int n, VD tau, VD K) {
    VD gradients = gradientVIXprice(p, r, n, tau, K, tbar);
    for (int j = 0; j < 30 * 5; j++) {
        std::cout << std::setprecision(16);
        std::cout << "gradient" << j << "   " << gradients[j] << std::endl;
    }
}

void showSPXcallPrices(modelPara mp, VD tarr, double S, VD karr, double r,
                       int n) {
    VD SPXprices = SPXprice(mp, tarr, S, karr, r, n);

    for (int j = 0; j < n; j++) {
        std::cout << std::fixed << std::setprecision(16);
        std::cout << "strike_price  " << karr[j] << std::endl;
        std::cout << "maturity      " << tarr[j] << std::endl;
        std::cout << "SPXcall_price " << SPXprices[j] << std::endl << std::endl;
    }
}

void printVIXcalls(modelPara p, VD tau, double tbar, VD K, double r, int n) {
    VD VIXcalls = VIXprice(p, tau, tbar, K, r, n);
    for (int j = 0; j < n; j++) {
        std::cout << std::setprecision(16);
        std::cout << "strike_price  " << K[j] << std::endl;
        std::cout << "maturity      " << tau[j] << std::endl;
        std::cout << "VIXcall_price " << VIXcalls[j] << std::endl << std::endl;
    }
}

void printCFvol(modelPara p, double tau, int n) {
    CFvolData CFv;
    VD u = *gl.abs;
    for (int j = 0; j < n; j++) {
        std::cout << std::setprecision(16);
        std::cout << "u " << u[j] << std::endl;
        CFv = CFvol(u[j] * one, p, tau);
        std::cout << "CFvol " << CFv.CFvol << std::endl;
    }
}

void printIntegrandVIXoption(modelPara p, double tau, double K, double tbar,
                             int n) {
    double VIXint;
    VD u = *gl.abs;
    for (int j = 0; j < n; j++) {
        std::cout << std::setprecision(16);
        //std::cout << "u " << u[j] << ", ";
        VIXint = VIXintegrand(u[j], p, tau, K, tbar).VIXint;
        std::cout << "VIXintegrand" << j << "  " << VIXint << std::endl;
    }
}
