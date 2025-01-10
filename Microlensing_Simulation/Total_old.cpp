#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <time.h>
#include <functional>
#include <algorithm>
#include <ctime>
#include <vector>
#include <thread>
#include "../spline.h"
#include <random>
#include"./Paper4GetPsi.h" //包含头文件
#include <iomanip>

using std::sin;
using std::cos;
using std::tan;
using std::exp; 
using std::pow;
using std::log;
using std::sqrt;
using std::atan;



using namespace std;

#define PI acos(-1)


//这个程序是将CE一年能观测到的所有宏观透镜一起分析，不区分2像、3像和4像，其中kappa star的选取也是分bin，根据20% fstar来算的。
int main()
{
    time_t TimeStart;
    TimeStart = time(NULL); //start time.
    //储存时间序列长度以及saddle点的x10和x20。
    vector<long> TimeLengthMinimum;
    int NumMinimum = 0;
    int NumSaddle = 0;
    vector<long> TimeLengthSaddle;
    vector<double> X10Saddle;
    vector<double> X20Saddle;
    char TimeLengthFileMinimum[100];
    char TimeLengthFileSaddle[100];
    char X10File[100];
    char X20File[100];
    double f_star = 0.8; //这里需要自己调，自己换。[0.2, 0.5, 0.8]
    double kappaStarSet[] = {0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 100};//{0.06, 0.06, 0.06};
    vector<double> kappaThreshold;
    for(int star_i = 0; star_i < sizeof(kappaStarSet)/sizeof(f_star) - 1; star_i ++ )
    {
        kappaThreshold.push_back((kappaStarSet[star_i] + kappaStarSet[star_i + 1]) / f_star / 2);
        cout << kappaThreshold[star_i] << endl;
    }


    sprintf(TimeLengthFileMinimum,"ResultMinimum/TimeLength%0.1f.bin",f_star);
    sprintf(TimeLengthFileSaddle,"ResultSaddle/TimeLength%0.1f.bin",f_star); 
    sprintf(X10File,"ResultSaddle/X10%0.1f.bin",f_star);
    sprintf(X20File,"ResultSaddle/X20%0.1f.bin",f_star);
    ofstream TL_Minimum;
    ofstream TL_Saddle;
    ofstream X10_Saddle;
    ofstream X20_Saddle;
    TL_Minimum.open(TimeLengthFileMinimum, std::ofstream::binary);
    TL_Saddle.open(TimeLengthFileSaddle, std::ofstream::binary);
    X10_Saddle.open(X10File, std::ofstream::binary);
    X20_Saddle.open(X20File, std::ofstream::binary);
    //微透镜场的配置
    int PrecisionFactor = 100;
    double MassPerSolarMass = 1;
    double M_sun = 2 * pow(10,30);
    double G = 6.67 * pow(10,-11);
    double c = 3 * pow(10,8);
    double FileImageResolution = 0.1;
    double ImageResolution = 0.01;
    long EspBoost = (double) FileImageResolution/ImageResolution;
    double MicroSkyLimit = 1000; //透镜平面大小
    long MicroFileLenRows = (2*MicroSkyLimit)/FileImageResolution; 
    cout << "MicroFileLenRows = " << MicroFileLenRows << endl;
    long MicroFileLenColumns = (2*MicroSkyLimit)/FileImageResolution;
    double y1, y2;
    y1 = 0;//0.4*sin(PI/8);//0;//
    y2 = 0;//0.4*cos(PI/8); //0;   //  source position
    //读取透镜的配置
    ifstream kappafile("/disk1/home/shanxk/work/Paper4_copy_CE/SampleResult/kappa.csv", ios::in);
    ifstream gammafile("/disk1/home/shanxk/work/Paper4_copy_CE/SampleResult/gamma.csv", ios::in);
    ifstream AcceptLensindexfile("/disk1/home/shanxk/work/Paper4_copy_CE/SampleResult/AcceptLensIndex1YearFeiXuMD14_multi37360.csv", ios::in);
    ifstream lenszfile("/disk1/home/shanxk/work/Paper4_copy_CE/SampleResult/lens_z.csv", ios::in);
    ifstream imagenumfile("/disk1/home/shanxk/work/Paper4_copy_CE/SampleResult/imagenumber.csv", ios::in); 
    ifstream SNROnlyMacrofile("/disk1/home/shanxk/work/Paper4_copy_CE/Lensed_SampleResult/SNR_network_only_macro.csv", ios::in);
    string linekappa;
    string linegamma;
    string lineindex;
    string linelensz;
    string lineimagenum;
    string linesnr;
    vector<double> lenskappa;
    vector<double> lensgamma;
    vector<double> lensindex;
    vector<double> lensz;
    vector<double> lensimagenum;
    vector<double> lenssnr;
    int index_num = 0;
    while (getline(kappafile, linekappa))
    {
        istringstream sin(linekappa);
        
        string kappa_tmp;
        while (getline(sin, kappa_tmp, ','))
        {
            lenskappa.push_back(stod(kappa_tmp));
        }
    }
    while (getline(gammafile, linegamma))
    {
        istringstream sin(linegamma);
        
        string gamma_tmp;
        while (getline(sin, gamma_tmp, ','))
        {
            lensgamma.push_back(stod(gamma_tmp));
        }
    }
    while (getline(AcceptLensindexfile, lineindex))
    {
        istringstream sin(lineindex);
        
        string index_tmp;
        while (getline(sin, index_tmp, ','))
        {
            lensindex.push_back(stod(index_tmp));
            index_num ++;
        }
    }
    while (getline(lenszfile, linelensz))
    {
        istringstream sin(linelensz);
        
        string lensz_tmp;
        while (getline(sin, lensz_tmp, ','))
        {
            lensz.push_back(stod(lensz_tmp));
        }
    }
    while (getline(imagenumfile, lineimagenum))
    {
        istringstream sin(lineimagenum);
        
        string imagenum_tmp;
        while (getline(sin, imagenum_tmp, ','))
        {
            lensimagenum.push_back(stod(imagenum_tmp));
        }
    }
    while (getline(SNROnlyMacrofile, linesnr))
    {
        istringstream sin(linesnr);
        
        string lenssnr_tmp;
        while (getline(sin, lenssnr_tmp, ','))
        {
            lenssnr.push_back(stod(lenssnr_tmp));
        }
    }
    
    //以上
    
    
    int DisSet[] = {0};//{3555, 4555, 5555};
    int DisSetX[] = {0};//{5277};//{2, 2, 2};//{3378, 4378, 5378};//
    int DisSetY[] = {0};//{3831};//{4132, 5132, 6132};//
    int IndexkappaStarSet = 0; 
    int IndexSumimage = 0;
    for(int ImageIndex = 0; ImageIndex < index_num; ImageIndex ++)
    {
        int IndexInAccept = lensindex[ImageIndex];
        double LensRedshift = lensz[ImageIndex];
        double coeffi = 4 * G * MassPerSolarMass * M_sun * (1 + LensRedshift) / pow(c, 3);
        int SaveIndex = 0;
        for(int subImageIndex = 0; subImageIndex < lensimagenum[ImageIndex]; subImageIndex ++ )
        {
            if(lenssnr[IndexSumimage] > 12)
            {
                
            
                double kappa = lenskappa[IndexSumimage];//0.462;
                double gamma = lensgamma[IndexSumimage];//0.462;
                int kappaStarThresholdIndex = 0;
                double kappaStarThresholdSelected;
                for(int kappaStar_i = 0; kappaStar_i < sizeof(kappaStarSet)/sizeof(f_star) - 1; kappaStar_i ++ )
                {
                    if(kappaThreshold[kappaStar_i] < kappa)
                    {
                        kappaStarThresholdIndex += 1;
                    }
                }
                if(kappaStarThresholdIndex == sizeof(kappaStarSet)/sizeof(f_star) - 1)
                {
                    kappaStarThresholdSelected = kappaStarSet[kappaStarThresholdIndex - 1];
                }
                else
                {
                    kappaStarThresholdSelected = kappaStarSet[kappaStarThresholdIndex];
                }
                cout << "kappa star selected = " << kappaStarThresholdSelected << endl;
                cout << "f_* x kappa = " << kappa * f_star << endl;
                cout << "kappa = " << kappa << endl;
                double mur = 1 - kappa + gamma;
                double mut = 1 - kappa - gamma;
                char AreaName[100];
                char TimeName[100];
                
                //Minimum
                if((mur > 0)&&(mut > 0))
                {
                    sprintf(AreaName,"ResultMinimum/Total%d_%dReadAreaMinimum%0.1f.bin", IndexInAccept, SaveIndex, f_star);
                    sprintf(TimeName,"ResultMinimum/Total%d_%dReadTimeMinimum%0.1f.bin", IndexInAccept, SaveIndex, f_star); 
                    ofstream fp_area;
                    ofstream fp_time;
                    fp_area.open(AreaName, std::ofstream::binary);
                    fp_time.open(TimeName, std::ofstream::binary);
                    
                    
                    
            
                    double Axisa = sqrt(1/coeffi/mur);
                    double Axisb = sqrt(1/coeffi/mut);
                    
                    /*变量*/
                    double* ReceivePre = Preparation4CreatPhiKappaStar("Minimum", kappa, gamma, kappaStarThresholdSelected, MassPerSolarMass, LensRedshift, PrecisionFactor);
                    double SkyLimit = ReceivePre[0];
                    int NStar = ReceivePre[1];
                    double TimeNeed2CalMax = ReceivePre[2];
                    delete [] ReceivePre;
                    /*天区坐标*/
                    long FileLengthX = (2*SkyLimit + FileImageResolution)/FileImageResolution;
                    cout << "Sky need to calculate length = " << FileLengthX << endl;
                    double* FileX1Set = new double [FileLengthX];
                    long FileX1Index = 0;
                    for(double TmpX1 = -SkyLimit + FileImageResolution/2 ; TmpX1 <= SkyLimit; TmpX1 += FileImageResolution)
                    {
                        if(FileX1Index < FileLengthX)
                        {
                            FileX1Set[FileX1Index] = TmpX1;
                            FileX1Index ++ ;
                        }
                        else
                        {
                            FileX1Index ++ ;
                        }
                    } 

                    if(FileX1Index == FileLengthX)
                    {
                        cout << "File X1's length is right "<< endl;
                        cout << "Length of FileX1 set = " << FileLengthX << " =  Maximum FileFX1 set index + 1 = " << FileX1Index << endl; 
                    }
                    else if(FileX1Index > FileLengthX)
                    {
                        cout << "Length of FileX1 set = " << FileLengthX << " <  Maximum FileFX1 set index + 1 = " << FileX1Index << endl; 
                        double FileX1Start = - SkyLimit + FileImageResolution/2;
                        for(long FileX1IndexTmp = 0 ; FileX1IndexTmp < FileLengthX; FileX1IndexTmp ++ )
                        {
                            FileX1Set[FileX1IndexTmp] = FileX1Start + FileX1IndexTmp * FileImageResolution;
                        
                        } 
                        FileX1Index = FileLengthX; //FIXME这里是我新加的，在运行0.4、0.6、0.8的时候没有这句话，但是我觉得这句话应该加上，好处是在上面三种情况运行的时候，没有一个事例运行过这部分，所以后面的Saddle部分，我没有修改。
                    }
                    else
                    {
                        cout << "Length of FileX1 set = " << FileLengthX << " >  Maximum FileFX1 set index + 1 = " << FileX1Index << endl; 
                        double FileX1Start = - SkyLimit + FileImageResolution/2;
                        for(long FileX1IndexTmp = 0 ; FileX1IndexTmp < FileLengthX; FileX1IndexTmp ++ )
                        {
                            FileX1Set[FileX1IndexTmp] = FileX1Start + FileX1IndexTmp * FileImageResolution;
                        
                        } 
                        FileX1Index = FileLengthX;

                    }    
                    long LengthX = 0;
                    vector<double> X1Set;
                    for(double TmpX1 = FileX1Set[0]; TmpX1 <= FileX1Set[FileX1Index - 1]; TmpX1 += ImageResolution)
                    {
                        X1Set.push_back(TmpX1);
                        LengthX += 1;
                    }
                    cout << "Boost Length X = " << LengthX << endl;
                    double* X2Set = new double [LengthX];
                    for(long TmpX2Index = 0; TmpX2Index < LengthX; TmpX2Index ++ )
                    {
                        X2Set[TmpX2Index] = X1Set[TmpX2Index];
                    }
                    cout << "X1Set[0] = " << X1Set[0] << endl;
                    cout << "X2Set[0] = " << X2Set[0] << endl; 
                    cout << "X2Set[-1] = " << X2Set[LengthX - 1] << endl;
                    cout << "X1Set[-1] = " << X1Set[LengthX - 1] << endl;
                    cout << "FileX1Set[0] = " << FileX1Set[0] << endl; 
                    cout << "FileX1Set[FileLengthX - 1] = " << FileX1Set[FileLengthX - 1] << endl;
                    cout << "FileX1Set[FileX1Index - 1] = " << FileX1Set[FileX1Index - 1] << endl;

                    
                    
                    

                    //判断有多大的范围可以移动
                    long DisColumnsMax = MicroFileLenColumns - FileLengthX; 
                    long DisRowsMax = MicroFileLenRows - FileLengthX;
                    cout << "DisColumnsmax = " << DisColumnsMax << endl;
                    cout << "DisRowsMax = " << DisRowsMax << endl;
                    if(DisColumnsMax != 0 and DisRowsMax != 0)
                    {
                        DisColumnsMax = rand()%DisColumnsMax;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = rand()%DisRowsMax;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else if(DisColumnsMax != 0 and DisRowsMax == 0)
                    {
                        DisColumnsMax = rand()%DisColumnsMax;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = 0;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else if(DisColumnsMax == 0 and DisRowsMax != 0)
                    {
                        DisColumnsMax = 0;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = rand()%DisRowsMax;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else
                    {
                        DisColumnsMax = 0;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = 0;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    
                    cout << "Random DisColumns = " << DisColumnsMax << endl;
                    cout << "Random DisRows = " << DisRowsMax << endl;
                    
                    double RudeImageResolution; // SkyLimit, x2: lens plane coordinate. 
                    //y1, y2: source plane coordinate. esp: lens plane resolution

                    long RudeEspBoost = EspBoost*10;
                    RudeImageResolution = ImageResolution * RudeEspBoost;

                    
                    
                    /*下面是粗略计算最大值点和最小值点的子程序*/
                    double* RudeTimeDelayMinAndMax = new double[2]; //这个是用来存后面得到的最大值和最小值的。
                    RudeTimeDelayMinAndMax[0] = 10000;
                    RudeTimeDelayMinAndMax[1] = -10000;
                    double TestX1, TestX2;
                    long TestX1Index, TestX2Index;
                    for(long RudeX1Index = 0; RudeX1Index < LengthX; RudeX1Index += RudeEspBoost )
                    {
                        double RudeX1Tmp = X1Set[RudeX1Index];

                        double* RudePsi = ReadCreatPsi("Minimum", kappa, gamma, kappaStarThresholdSelected, LengthX, MicroFileLenRows, MicroFileLenColumns, FileLengthX, FileX1Set, X2Set, RudeX1Index, DisColumnsMax, DisRowsMax, y1, y2, EspBoost);
                        // cout << RudeX1Tmp << endl;
                        for(long RudeY1Index = 0; RudeY1Index < LengthX; RudeY1Index += RudeEspBoost )
                        {
                            double RudeX2Tmp = X2Set[RudeY1Index];
                            double MacroT_tmp = pow(RudeX1Tmp - y1/mur,2)/(2*Axisa*Axisa) + pow(RudeX2Tmp - y2/mut,2)/(2*Axisb*Axisb);
                            if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                            {
                                if(RudePsi[RudeY1Index] > RudeTimeDelayMinAndMax[1])
                                {
                                    RudeTimeDelayMinAndMax[1] = RudePsi[RudeY1Index];
                                    
                                }
                                if(RudePsi[RudeY1Index] < RudeTimeDelayMinAndMax[0])
                                {
                                    RudeTimeDelayMinAndMax[0] = RudePsi[RudeY1Index];
                                    TestX1 = X1Set[RudeX1Index];
                                    TestX2 = X2Set[RudeY1Index];
                                    TestX1Index = RudeX1Index;
                                    TestX2Index = RudeY1Index;
                                    // cout << "Test = " << X1Set[TestX1Index] << endl;
                                    
                                }
                            }
                        }
                        delete [] RudePsi;
                    }
                    double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0] * coeffi;
                    double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1] * coeffi;
                    cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
                    cout << "Rude Minimum X1 = " << TestX1 << " X2 = " << TestX2 << endl;
                    cout << "Rude Minimum X1 Index = " << TestX1Index << " X2 Index = " << TestX2Index << endl; 
                    cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
                    delete[] RudeTimeDelayMinAndMax;
                    /*以上*/

                    /*下面并行计算时间延迟曲线*/

                    //the next is for time range
                    double TimeDelayStart = RudeTimeDelayMin; 
                    double TimeDelayEsp = pow(10,-6);
                    double TimeDelayEnd = RudeTimeDelayMax;
                    //
                    cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
                    long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
                    // the length of time delay area
                    cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
                    TimeLengthMinimum.push_back(LenTimeDelay);
                    NumMinimum += 1;
                    
                    double* TimeDelayRange = new double [LenTimeDelay];
                    double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
                    //generate time array.
                    for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
                    {
                    
                        TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
                        TimeDelayTmp += TimeDelayEsp;
                        
                        

                    }
                    cout << TimeDelayRange[LenTimeDelay - 1] << endl;
                    //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


                    long* TimeDelayArea = new long [LenTimeDelay](); //存放时间延迟面积的数组。


                    int ThreadCount = 100; //使用的线程数量。
                    std::thread* threads[ThreadCount];
                    int args[ThreadCount];
                    //generate SkyLimit array
                    
                    /*
                    double* X1Set = new double [LenX1Set];
                    double X1Tmp = -SkyLimit;
                    for(long X1Index = 0; X1Index < LenX1Set; X1Index ++ ) //创建像平面的坐标。
                    {
                    
                        X1Set[X1Index] = X1Tmp;
                        X1Tmp += RudeImageResolution;
                        
                    }
                    //cout << "X1Set[-1]:" << X1Set[LenX1Set - 1] << endl;
                    */

                    /*用来存放低于粗估计的时间延迟最小值或者高于粗估计的时间延迟最大值的数据*/
                    double* TimeDelayBelow = new double [ThreadCount]();
                    long TimeDelayBelowIndex = 0;
                    double* TimeDelayUpper = new double [ThreadCount]();
                    long TimeDelayUpperIndex = 0;
                    
                    //创建二维的数组用来存放时间延迟面积的数据，并初始化为0；
                    long** TimeDelayArea2D = new long*[ThreadCount];
                    for(int i = 0; i < ThreadCount; i++)
                    {
                        TimeDelayArea2D[i] = new long[LenTimeDelay]();
                    }
                    //存放cut以后的微透镜加负质量片场
                    
                    double** MicroAndMinus2D = new double*[LengthX];
                    for(long i = 0; i < LengthX; i ++ )
                    {
                        MicroAndMinus2D[i] = new double[LengthX];
                    }
                
                    

                    //倍乘因子，主要是算一下一个核算多少行。
                    long Order = LengthX/ThreadCount; //这里选取线程数的时候，最好可以整除，要不然会变慢很多。
                    cout << "One Core Has " << Order << " Task" << endl;

                    auto task = [&] (int i) -> void {

                        vector<double> TimeDelayPrecisionMin; //存放实验延迟最小的精确解。
                        vector<double> TimeDelayPrecisionMax;
                        long CountTimeDelayLower = 0; //数一下有多少时间延迟小于了粗略解。
                        long CountTimeDelayUpper = 0;
                    
                        long IndexUpLimit = i + 1 == ThreadCount ? LengthX : (i + 1) * Order; //三目运算符
                        for(long Index = i*Order; Index < IndexUpLimit; Index++ )
                        { 
                            double TmpX1 = X1Set[Index]; //算一下这一列的x值。
                            double* TimeDelay = ReadCreatPsi("Minimum", kappa, gamma, kappaStarThresholdSelected, LengthX, MicroFileLenRows, MicroFileLenColumns, FileLengthX, FileX1Set, X2Set, Index, DisColumnsMax, DisRowsMax, y1, y2, EspBoost);
                            
                            
                            // double* TimeDelay = CreatPsi(kappa, gamma, TmpX1, MicroLensCoorXY, NStar, SkyLimit, LengthX, X1Set, kappaStar, y1, y2);
                            for(long JIndex = 0; JIndex < LengthX; JIndex += 1)
                            {
                                double TmpX2 = X2Set[JIndex]; //算一下坐标的y值。
                                double MacroT_tmp = pow(TmpX1 - y1/mur,2)/(2*Axisa*Axisa) + pow(TmpX2 - y2/mut,2)/(2*Axisb*Axisb);
                                if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                                //判断是不是坐标是不是在双曲线区域内。
                                // if(1)
                                {
                                    double TimeDelayTmp = TimeDelay[JIndex] * coeffi; // add time delay result to result array.
                                    
                                    if(TimeDelayTmp < RudeTimeDelayMin)
                                    {
                                        TimeDelayPrecisionMin.push_back(TimeDelayTmp);
                                        CountTimeDelayLower += 1;
                                    }
                                    else if(TimeDelayTmp > RudeTimeDelayMax)
                                    {
                                        TimeDelayPrecisionMax.push_back(TimeDelayTmp);
                                        CountTimeDelayUpper += 1;
                                    }
                                    else
                                    {
                                        long TimeDelayIndex = long ((TimeDelayTmp - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 1;
                                        // cout << time_index << endl;
                                        
                                        TimeDelayArea2D[i][TimeDelayIndex] += 1; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                    
                                    }
                                }

                            }
                            delete[] TimeDelay;
                        
                        }
                        
                        if(CountTimeDelayLower != 0)
                        {
                        double TimeDelayPrecisionMinMin = *min_element(TimeDelayPrecisionMin.begin(), TimeDelayPrecisionMin.end()); 
                        TimeDelayBelow[i] = TimeDelayPrecisionMinMin;
                        }
                        else
                        {
                            TimeDelayBelow[i] = RudeTimeDelayMin;
                        }

                        if(CountTimeDelayUpper != 0)
                        {
                        double TimeDelayPrecisionMaxMax = *max_element(TimeDelayPrecisionMax.begin(), TimeDelayPrecisionMax.end());
                        TimeDelayUpper[i] = TimeDelayPrecisionMaxMax; 
                        }
                        else
                        {
                            TimeDelayUpper[i] = 0;
                        }
                        

                    };

                    auto run_tasks = [&] (int count) -> void {
                        for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                            threads[CountIndex] = new std::thread (task, args[CountIndex]);
                        }
                        for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                            threads[CountIndex]->join();
                            delete threads[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
                        }
                    };

                    long JThread = 0;
                    for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
                    {
                        args[JThread++] = ThreadIndex;
                        if (JThread == ThreadCount) {
                            run_tasks(ThreadCount);
                            JThread = 0;
                        }
                    }
                    run_tasks(JThread);
                    
                    for(long AreaIndex = 0; AreaIndex < LenTimeDelay; AreaIndex ++)
                    {
                        for(long AreaSumIndex = 0; AreaSumIndex < ThreadCount; AreaSumIndex ++)
                        {
                            TimeDelayArea[AreaIndex] += TimeDelayArea2D[AreaSumIndex][AreaIndex];
                        }
                    }

                    
                    

                    //根据最小的时间，重标度一下时间。
                    double DeltaFinalPrecisionTime = RudeTimeDelayMin - *min_element(TimeDelayBelow, TimeDelayBelow + ThreadCount);
                    
                    cout << "Delta Time = " << DeltaFinalPrecisionTime << endl;
                    if(DeltaFinalPrecisionTime < 0)
                    {
                        cout << "Don't include the minimum point in this Hyperbola area" << endl;
                        DeltaFinalPrecisionTime = 0;
                    }
                    fp_area.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(long)*LenTimeDelay);
                    fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
                    
                    fp_area.close();
                    fp_time.close();
                    
                    cout << "Precision Minimum Time Delay =" << *min_element(TimeDelayBelow, TimeDelayBelow + ThreadCount) << endl;
                    cout << "Precision Maximum Time Delay =" << *max_element(TimeDelayUpper, TimeDelayUpper + ThreadCount) << endl;

                    //print run time.
                    double CostTime = time(NULL) - TimeStart;
                    cout << "all time : " << CostTime  << " s" << endl;

                    delete []TimeDelayRange;
                    delete []TimeDelayArea;
                    delete []TimeDelayArea2D;
                    delete []TimeDelayBelow;
                    delete []TimeDelayUpper;
                    //delete []X1Set;
                    //delete[] TimeDelay;
                }

                //Saddle
                
                if(mur * mut < 0)
                {
                    mur = 1 - kappa + gamma;
                    mut = kappa + gamma - 1; //FIXME注意，这里有可能反号，好处是在SIE模型中kappa = gamma，所以mur = 1, mut<0
                    if(mur < 0)
                    {
                        cout << "*********************************************" << endl;
                        cout << "WRONG Type II kappa and shear !!!!!!!!!!!!!!!" << endl;
                    }
                    sprintf(AreaName,"ResultSaddle/Total%d_%dReadAreaSaddle%0.1f.bin", IndexInAccept, SaveIndex, f_star);
                    sprintf(TimeName,"ResultSaddle/Total%d_%dReadTimeSaddle%0.1f.bin", IndexInAccept, SaveIndex, f_star);  
                    ofstream fp_area;
                    ofstream fp_time;
                    fp_area.open(AreaName, std::ofstream::binary);
                    fp_time.open(TimeName, std::ofstream::binary);
                    
                    double Axisa = sqrt(1/coeffi/mur);
                    double Axisb = sqrt(1/coeffi/mut);
                    
                    /*变量*/
                    double* ReceivePre = SaddlePreparation4CreatPhiKappaStar("Saddle", kappa, gamma, kappaStarThresholdSelected, MassPerSolarMass, LensRedshift, PrecisionFactor);
                    double X10New = ReceivePre[0];
                    double X20New = ReceivePre[6];
                    double Epsilon1 = ReceivePre[7];
                    double Epsilon2 = ReceivePre[8];
                    int NStar = ReceivePre[1];
                    double TimeNeed2CalMax = ReceivePre[2];
                    double SkyLimitX = ReceivePre[3];
                    double SkyLimitY = ReceivePre[4];
                    double kappaStarNew = ReceivePre[5];
                    cout << "kappaStarNew = " << kappaStarNew << endl;
                    cout << "X10New = " << X10New << "X20New = " << X20New << endl;
                    X10Saddle.push_back(X10New);
                    X20Saddle.push_back(X20New);
                    
                    delete [] ReceivePre;
                    long FileLengthX = (2*SkyLimitX + FileImageResolution)/FileImageResolution;
                    long FileLengthY = (2*SkyLimitY + FileImageResolution)/FileImageResolution;
                    /*天区横坐标*/
                    /*文件中的*/
                    double* FileX1Set = new double [FileLengthX];
                    long FileX1Index = 0;
                    for(double TmpX1 = -SkyLimitX + FileImageResolution/2 ; TmpX1 <= SkyLimitX; TmpX1 += FileImageResolution)
                    {
                        if(FileX1Index < FileLengthX)
                        {
                            FileX1Set[FileX1Index] = TmpX1;
                            FileX1Index ++ ;
                        }
                        else
                        {
                            FileX1Index ++ ;
                        }
                    } 
                    if(FileX1Index == FileLengthX)
                    {
                        cout << "FileX1's length is right "<< endl;
                        cout << "Length of FileX1 set = " << FileLengthX << " =  Maximum FileX1 set index + 1 = " << FileX1Index << endl; 
                    }
                    else if(FileX1Index > FileLengthX)
                    {
                        cout << "Length of FileX1 set = " << FileLengthX << " <  Maximum FileFX1 set index + 1 = " << FileX1Index << endl; 
                        double FileX1Start = - SkyLimitX + FileImageResolution/2;
                        for(long FileX1IndexTmp = 0 ; FileX1IndexTmp < FileLengthX; FileX1IndexTmp ++ )
                        {
                            FileX1Set[FileX1IndexTmp] = FileX1Start + FileX1IndexTmp * FileImageResolution;
                        
                        } 
                    }
                    else
                    {
                        cout << "Length of FileX1 set = " << FileLengthX << " >  Maximum FileFX1 set index + 1 = " << FileX1Index << endl; 
                        double FileX1Start = - SkyLimitX + FileImageResolution/2;
                        for(long FileX1IndexTmp = 0 ; FileX1IndexTmp < FileLengthX; FileX1IndexTmp ++ )
                        {
                            FileX1Set[FileX1IndexTmp] = FileX1Start + FileX1IndexTmp * FileImageResolution;
                        
                        } 
                        FileX1Index = FileLengthX;

                    }   
                    /*新生成的*/
                    long LengthX = 0;
                    vector <double> X1SetVec;
                    for(double TmpX1 = FileX1Set[0]; TmpX1 <= FileX1Set[FileX1Index - 1]; TmpX1 += ImageResolution)
                    {
                        X1SetVec.push_back(TmpX1);
                        LengthX += 1;
                    }
                    double * X1Set = new double [LengthX];
                    for(long TmpX1Index = 0; TmpX1Index < LengthX; TmpX1Index ++ )
                    {
                        X1Set[TmpX1Index] = X1SetVec[TmpX1Index];
                    }
                    
                    cout << "X1Set[-1] = " << X1Set[LengthX - 1] << endl;
                    cout << "FileX1Set[-1] = " << FileX1Set[FileLengthX - 1] << endl;
                    
                    
                    
                    //天区纵坐标
                    /*文件中的*/
                    double* FileX2Set = new double [FileLengthY];
                    long FileX2Index = 0;
                    for(double TmpX2 = -SkyLimitY + FileImageResolution/2 ; TmpX2 <= SkyLimitY; TmpX2 += FileImageResolution)
                    {
                        if(FileX2Index < FileLengthY)
                        {
                            FileX2Set[FileX2Index] = TmpX2;
                            FileX2Index ++ ;
                        }
                        else
                        {
                            FileX2Index ++ ;
                        }
                    } 
                    if(FileX2Index == FileLengthY)
                    {
                        cout << "FileX2's length is right "<< endl;
                        cout << "Length of FileX2 set = " << FileLengthY << " =  Maximum FileX2 set index + 1 = " << FileX2Index << endl; 
                    }
                    else if(FileX2Index > FileLengthY)
                    {
                        cout << "Length of FileX2 set = " << FileLengthY << " <  Maximum FileX2 set index + 1 = " << FileX2Index << endl; 
                        double FileX2Start = - SkyLimitY + FileImageResolution/2;
                        for(long FileX2IndexTmp = 0 ; FileX2IndexTmp < FileLengthY; FileX2IndexTmp ++ )
                        {
                            FileX2Set[FileX2IndexTmp] = FileX2Start + FileX2IndexTmp * FileImageResolution;
                        
                        } 
                    }
                    else
                    {
                        cout << "Length of FileX2 set = " << FileLengthY << " >  Maximum FileX2 set index + 1 = " << FileX2Index << endl; 
                        double FileX2Start = - SkyLimitY + FileImageResolution/2;
                        for(long FileX2IndexTmp = 0 ; FileX2IndexTmp < FileLengthY; FileX2IndexTmp ++ )
                        {
                            FileX2Set[FileX2IndexTmp] = FileX2Start + FileX2IndexTmp * FileImageResolution;
                        
                        } 
                        FileX2Index = FileLengthY;

                    }   
                    /*新生成的*/
                    long LengthY = 0;
                    vector <double> X2SetVec;
                    for(double TmpX2 = FileX2Set[0]; TmpX2 <= FileX2Set[FileX2Index - 1]; TmpX2 += ImageResolution)
                    {
                        X2SetVec.push_back(TmpX2);
                        LengthY += 1;
                    }
                    double * X2Set = new double [LengthY];
                    for(long TmpX2Index = 0; TmpX2Index < LengthY; TmpX2Index ++ )
                    {
                        X2Set[TmpX2Index] = X2SetVec[TmpX2Index];
                    }
                    
                    cout << "X2Set[-1] = " << X2Set[LengthY - 1] << endl;
                    cout << "FileX2Set[-1] = " << FileX2Set[FileLengthY - 1] << endl;
                    
                    cout << "LengthX = " << LengthX << endl;
                    cout << "LengthY = " << LengthY << endl;

                    //判断有多大的范围可以移动
                    long DisColumnsMax = MicroFileLenColumns - FileLengthX; 
                    long DisRowsMax = MicroFileLenRows - FileLengthY;
                    cout << "DisColumnsmax = " << DisColumnsMax << endl;
                    cout << "DisRowsMax = " << DisRowsMax << endl;
                    if(DisColumnsMax != 0 and DisRowsMax != 0)
                    {
                        DisColumnsMax = rand()%DisColumnsMax;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = rand()%DisRowsMax;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else if(DisColumnsMax != 0 and DisRowsMax == 0)
                    {
                        DisColumnsMax = rand()%DisColumnsMax;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = 0;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else if(DisColumnsMax == 0 and DisRowsMax != 0)
                    {
                        DisColumnsMax = 0;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = rand()%DisRowsMax;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    else
                    {
                        DisColumnsMax = 0;// DisColumnsMax / 2;//DisSet[IndexkappaStarSet];// DisColumnsMax/DisSet[IndexkappaStarSet];//0;//
                        DisRowsMax = 0;//DisRowsMax / 2;//DisSet[IndexkappaStarSet];// DisRowsMax/DisSet[IndexkappaStarSet];//0;//
                    }
                    cout << "Random DisColumns = " << DisColumnsMax << endl;
                    cout << "Random DisRows = " << DisRowsMax << endl;

                    
                    double RudeImageResolution; // SkyLimit, x2: lens plane coordinate. 
                    //y1, y2: source plane coordinate. esp: lens plane resolution.

                    long RudeEspBoost = EspBoost*10;
                    cout << "RudeEspBoost = " << RudeEspBoost << endl;
                    RudeImageResolution = ImageResolution * RudeEspBoost;
                    
                    /*下面是粗略计算最大值点和最小值点的子程序*/
                    double* RudeTimeDelayMinAndMax = new double[2]; //这个是用来存后面得到的最大值和最小值的。
                    RudeTimeDelayMinAndMax[0] = 1000;
                    RudeTimeDelayMinAndMax[1] = -1000;
                    double TestX1, TestX2;
                    long TestX1Index, TestX2Index;
                    for(long RudeX1Index = 0; RudeX1Index < LengthX; RudeX1Index += RudeEspBoost )
                    {

                        double RudeX1Tmp = X1Set[RudeX1Index];
                        double* RudePsi = SaddleReadCreatPsi(kappa, gamma, kappaStarThresholdSelected, LengthY, MicroFileLenRows, MicroFileLenColumns, FileLengthY, FileX2Set, X1Set, X2Set, RudeX1Index, DisColumnsMax, DisRowsMax, y1, y2, EspBoost);
                        
                        
                        // cout << RudeX1Tmp << endl;
                        for(long RudeY1Index = 0; RudeY1Index < LengthY; RudeY1Index += RudeEspBoost )
                        {
                            double RudeX2Tmp = X2Set[RudeY1Index];
                            double MacroT_tmp = pow(RudeX1Tmp - y1/mur,2)/(2*Axisa*Axisa) - pow(RudeX2Tmp - y2/mut,2)/(2*Axisb*Axisb);
                            if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                            {
                                if(MacroT_tmp >= 0)
                                {
                                    if((RudeX2Tmp > - X20New)&&(RudeX2Tmp < X20New))
                                    // if(1)
                                    {
                                        double RudeTimeDelayTmp = RudePsi[RudeY1Index] * coeffi;
                                        if(RudeTimeDelayTmp > RudeTimeDelayMinAndMax[1])
                                        {
                                            RudeTimeDelayMinAndMax[1] = RudeTimeDelayTmp;
                                        }
                                        if(RudeTimeDelayTmp < RudeTimeDelayMinAndMax[0])
                                        {
                                            RudeTimeDelayMinAndMax[0] = RudeTimeDelayTmp;
                                        }
                                    }
                                }

                                else if(MacroT_tmp < 0)
                                {
                                    if((RudeX1Tmp > - X10New)&&(RudeX1Tmp < X10New))
                                    {
                                        double RudeTimeDelayTmp = RudePsi[RudeY1Index] * coeffi;
                                        if(RudeTimeDelayTmp > RudeTimeDelayMinAndMax[1])
                                        {
                                            RudeTimeDelayMinAndMax[1] = RudeTimeDelayTmp;
                                        }
                                        if(RudeTimeDelayTmp < RudeTimeDelayMinAndMax[0])
                                        {
                                            RudeTimeDelayMinAndMax[0] = RudeTimeDelayTmp;
                                        }
                                    }
                                }
                            }
                        }
                        delete [] RudePsi;
                    }



                    double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0];
                    double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1];
                    cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
                    cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
                    delete[] RudeTimeDelayMinAndMax;
                    /*以上*/

                    /*下面并行计算时间延迟曲线*/

                    //the next is for time range
                    double TimeDelayStart = RudeTimeDelayMin; 
                    double TimeDelayEsp = pow(10,-6);
                    double TimeDelayEnd = RudeTimeDelayMax;
                    //
                    cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
                    long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
                    // the length of time delay area
                    cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
                    TimeLengthSaddle.push_back(LenTimeDelay);
                    NumSaddle += 1;
                    
                    double* TimeDelayRange = new double [LenTimeDelay];
                    double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
                    //generate time array.
                    for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
                    {
                    
                        TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
                        TimeDelayTmp += TimeDelayEsp;
                        
                        

                    }
                    cout << TimeDelayRange[LenTimeDelay - 1] << endl;
                    //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


                    long* TimeDelayArea = new long [LenTimeDelay](); //存放时间延迟面积的数组。


                    int ThreadCount = 100; //使用的线程数量。
                    std::thread* threads[ThreadCount];
                    int args[ThreadCount];
                    //generate SkyLimit array
                    
                    /*
                    double* X1Set = new double [LenX1Set];
                    double X1Tmp = -SkyLimit;
                    for(long X1Index = 0; X1Index < LenX1Set; X1Index ++ ) //创建像平面的坐标。
                    {
                    
                        X1Set[X1Index] = X1Tmp;
                        X1Tmp += RudeImageResolution;
                        
                    }
                    //cout << "X1Set[-1]:" << X1Set[LenX1Set - 1] << endl;
                    */

                    /*用来存放低于粗估计的时间延迟最小值或者高于粗估计的时间延迟最大值的数据*/
                    double* TimeDelayBelow = new double [ThreadCount]();
                    long TimeDelayBelowIndex = 0;
                    double* TimeDelayUpper = new double [ThreadCount]();
                    long TimeDelayUpperIndex = 0;
                    
                    //创建二维的数组用来存放时间延迟面积的数据，并初始化为0；
                    long** TimeDelayArea2D = new long*[ThreadCount];
                    for(int i = 0; i < ThreadCount; i++)
                    {
                        TimeDelayArea2D[i] = new long[LenTimeDelay]();
                    }
                    //存放cut以后的微透镜加负质量片场
                    
                    double** MicroAndMinus2D = new double*[LengthX];
                    for(long i = 0; i < LengthX; i ++ )
                    {
                        MicroAndMinus2D[i] = new double[LengthY];
                    }

                    //倍乘因子，主要是算一下一个核算多少行。
                    long Order = LengthX/ThreadCount; //这里选取线程数的时候，最好可以整除，要不然会变慢很多。
                    cout << "One Core Has " << Order << " Task" << endl;

                    auto task = [&] (int i) -> void {

                        vector<double> TimeDelayPrecisionMin; //存放实验延迟最小的精确解。
                        vector<double> TimeDelayPrecisionMax;
                        long CountTimeDelayLower = 0; //数一下有多少时间延迟小于了粗略解。
                        long CountTimeDelayUpper = 0;
                    
                        long IndexUpLimit = i + 1 == ThreadCount ? LengthX : (i + 1) * Order; //三目运算符
                        for(long Index = i*Order; Index < IndexUpLimit; Index++ )
                        {
                            double TmpX1 = X1Set[Index]; //算一下这一列的x值。
                            double* TimeDelay = SaddleReadCreatPsi(kappa, gamma, kappaStarThresholdSelected, LengthY, MicroFileLenRows, MicroFileLenColumns, FileLengthY, FileX2Set, X1Set, X2Set, Index, DisColumnsMax, DisRowsMax, y1, y2, EspBoost);
                            
                        
                            // double* TimeDelay = CreatPsi(kappa, gamma, TmpX1, MicroLensCoorXY, NStar, SkyLimitX, SkyLimitY, LengthX, LengthY, X1Set, X2Set, kappaStar, y1, y2);
                            // CreatPsi(kappa, gamma, TmpX1, MicroLensCoorXY, NStar, SkyLimit, LengthX, X1Set, kappaStar);
                            for(long JIndex = 0; JIndex < LengthY; JIndex += 1)
                            {
                                double TmpX2 = X2Set[JIndex]; //算一下坐标的y值。
                                double MacroT_tmp = pow(TmpX1 - y1/mur,2)/(2*Axisa*Axisa) - pow(TmpX2 + y2/mut,2)/(2*Axisb*Axisb);
                                if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                                //判断是不是坐标是不是在双曲线区域内。
                                // if(1)
                                {
                                    if(MacroT_tmp >=0)
                                    {
                                        if((TmpX2 > - X20New)&&(TmpX2 < X20New))
                                        {
                                            double TimeDelayTmp = TimeDelay[JIndex] * coeffi; // add time delay result to result array.
                                            
                                            if(TimeDelayTmp <= RudeTimeDelayMin)
                                            {
                                                TimeDelayPrecisionMin.push_back(TimeDelayTmp);
                                                CountTimeDelayLower += 1;
                                            }
                                            else if(TimeDelayTmp >= RudeTimeDelayMax)
                                            {
                                                TimeDelayPrecisionMax.push_back(TimeDelayTmp);
                                                CountTimeDelayUpper += 1;
                                            }
                                            else
                                            {
                                                long TimeDelayIndex = long ((TimeDelayTmp - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 1;
                                                // cout << time_index << endl;
                                                if(TimeDelayIndex < 0)
                                                {
                                                    cout << "Error TimeDelayIndex = " << TimeDelayIndex << endl;
                                                }
                                                
                                                TimeDelayArea2D[i][TimeDelayIndex] += 1; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                            
                                            }
                                        }
                                    }
                                    else if(MacroT_tmp < 0)
                                    {
                                        if((TmpX1 > - X10New)&&(TmpX1 < X10New))
                                        {
                                            double TimeDelayTmp = TimeDelay[JIndex] * coeffi; // add time delay result to result array.
                                            
                                            if(TimeDelayTmp < RudeTimeDelayMin)
                                            {
                                                TimeDelayPrecisionMin.push_back(TimeDelayTmp);
                                                CountTimeDelayLower += 1;
                                            }
                                            else if(TimeDelayTmp > RudeTimeDelayMax)
                                            {
                                                TimeDelayPrecisionMax.push_back(TimeDelayTmp);
                                                CountTimeDelayUpper += 1;
                                            }
                                            else
                                            {
                                                long TimeDelayIndex = long ((TimeDelayTmp - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 1;
                                                // cout << time_index << endl;
                                                
                                                TimeDelayArea2D[i][TimeDelayIndex] += 1; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                            
                                            }
                                        }
                                    }
                                }

                            }
                            delete[] TimeDelay;
                        
                        }
                        
                        if(CountTimeDelayLower != 0)
                        {
                        double TimeDelayPrecisionMinMin = *min_element(TimeDelayPrecisionMin.begin(), TimeDelayPrecisionMin.end()); 
                        TimeDelayBelow[i] = TimeDelayPrecisionMinMin;
                        }
                        else
                        {
                            TimeDelayBelow[i] = RudeTimeDelayMin;
                        }

                        if(CountTimeDelayUpper != 0)
                        {
                        double TimeDelayPrecisionMaxMax = *max_element(TimeDelayPrecisionMax.begin(), TimeDelayPrecisionMax.end());
                        TimeDelayUpper[i] = TimeDelayPrecisionMaxMax; 
                        }
                        else
                        {
                            TimeDelayUpper[i] = 0;
                        }
                        

                    };

                    auto run_tasks = [&] (int count) -> void {
                        for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                            threads[CountIndex] = new std::thread (task, args[CountIndex]);
                        }
                        for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                            threads[CountIndex]->join();
                            delete threads[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
                        }
                    };

                    long JThread = 0;
                    for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
                    {
                        args[JThread++] = ThreadIndex;
                        if (JThread == ThreadCount) {
                            run_tasks(ThreadCount);
                            JThread = 0;
                        }
                    }
                    run_tasks(JThread);
                    
                    for(long AreaIndex = 0; AreaIndex < LenTimeDelay; AreaIndex ++)
                    {
                        for(long AreaSumIndex = 0; AreaSumIndex < ThreadCount; AreaSumIndex ++)
                        {
                            TimeDelayArea[AreaIndex] += TimeDelayArea2D[AreaSumIndex][AreaIndex];
                        }
                    }
                    
                    //根据最小的时间，重标度一下时间。
                    double DeltaFinalPrecisionTime = RudeTimeDelayMin - *min_element(TimeDelayBelow, TimeDelayBelow + ThreadCount);
                    
                    cout << "Delta Time = " << DeltaFinalPrecisionTime << endl;
                    if(DeltaFinalPrecisionTime < 0)
                    {
                        cout << "Don't include the minimum point in this Hyperbola area" << endl;
                        DeltaFinalPrecisionTime = 0;
                    }
                    fp_area.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(long)*LenTimeDelay);
                    fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
                    
                    fp_area.close();
                    fp_time.close();
                    
                    cout << "Precision Minimum Time Delay =" << *min_element(TimeDelayBelow, TimeDelayBelow + ThreadCount) << endl;
                    cout << "Precision Maximum Time Delay =" << *max_element(TimeDelayUpper, TimeDelayUpper + ThreadCount) << endl;

                    //print run time.
                    double CostTime = time(NULL) - TimeStart;
                    cout << "all time : " << CostTime  << " s" << endl;

                    delete []TimeDelayRange;
                    delete []TimeDelayArea;
                    delete []TimeDelayArea2D;
                    delete []TimeDelayBelow;
                    delete []TimeDelayUpper;
                    //delete []X1Set;
                    //delete[] TimeDelay;   
                }
                SaveIndex += 1;
            }
            IndexSumimage += 1;
        }
    }
    TL_Minimum.write(reinterpret_cast<char *>(&TimeLengthMinimum[0]), sizeof(long)*NumMinimum);
    TL_Saddle.write(reinterpret_cast<char *>(&TimeLengthSaddle[0]), sizeof(long)*NumSaddle);
    X10_Saddle.write(reinterpret_cast<char *>(&X10Saddle[0]), sizeof(double)*NumSaddle);
    X20_Saddle.write(reinterpret_cast<char *>(&X20Saddle[0]), sizeof(double)*NumSaddle);

    cout << "Number of Minimum" << NumMinimum << endl;
    cout << "Number of Saddle" << NumSaddle << endl;

    
    TL_Minimum.close();
    TL_Saddle.close();
    X10_Saddle.close();
    X20_Saddle.close();
    

    return 0;
}