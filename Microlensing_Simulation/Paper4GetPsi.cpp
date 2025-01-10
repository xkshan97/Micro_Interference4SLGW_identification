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
#include <string>
#include <random>
#include "../spline.h"

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

//This is for multi micro lenses, one need to set the \kappa_*.
double* Preparation4CreatPhiKappaStar(string MacroType, double kappa, double gamma, double kappaStar, double MassPerSolarMass, double LensRedshift, int PrecisionFactor)
{
    double* PreOutPut = new double[6]; //这个是用来存后面得到的最大值和最小值的。
    double MSun = 2 * pow(10,30);
    double GravityConstant = 6.67 * pow(10,-11);
    double LightSpeed = 3 * pow(10,8);
    double MicroLensCoeffi = 4 * GravityConstant * MassPerSolarMass * MSun * (1 + LensRedshift) / pow(LightSpeed, 3);
    double RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)});
    while(RNeed2Cal>900) //这里写的是900，为了留出来200的空隙
    {
        PrecisionFactor -= 1;
        RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}); 
    }
    // cout << "SNR = 60 exceed sky map" << endl;
    cout << "NOW SNR = " << PrecisionFactor << endl;
    double TimeNeed2CalMax = MicroLensCoeffi/2*min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)})*pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2);
    cout << "TimeNeed2CalMax = " << TimeNeed2CalMax << endl;
    double MagnificationR = 1 - kappa + gamma;
    double MagnificationT = 1 - kappa - gamma;
    double SkyLimit;
    double X20=0;
    double X10=0;
    if(MacroType == "Minimum")
    {
        double X10 = sqrt(2*TimeNeed2CalMax/MicroLensCoeffi/MagnificationR);
        double X20 = sqrt(2*TimeNeed2CalMax/MicroLensCoeffi/MagnificationT);
        SkyLimit = max({X10,X20}); 
        // if(SkyLimit<50)
        // {
        //     SkyLimit = 50;
        // }
        // cout << TmpX1 << TmpX2 << SkyLimit << endl;

    }
    
    
    else if(MacroType == "Maximum")
    {
        double X10 = sqrt(-2*TimeNeed2CalMax/MicroLensCoeffi/MagnificationR);
        double X20 = sqrt(-2*TimeNeed2CalMax/MicroLensCoeffi/MagnificationT);
        SkyLimit = max({X10,X20}); 
        // if(SkyLimit<50)
        // {
        //     SkyLimit = 50;
        // }
        // cout << TmpX1 << TmpX2 << SkyLimit << endl;

    }
    else
    {
        cout << "Unknown Macrolensing Type" << endl;
    }
    
    // if(SkyLimit < 100)
    // {
    //     SkyLimit = 100;
    // }
    
    double AreaS = 4 * pow(SkyLimit, 2); //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
    int NStar = kappaStar * AreaS / PI;
    double kappaStarNew = NStar * PI / AreaS;
    cout << "Nstar = " << NStar << endl;
    cout << "SkyLimit = " << SkyLimit << endl;
    PreOutPut[0] = SkyLimit;
    PreOutPut[1] = NStar;
    PreOutPut[2] = TimeNeed2CalMax;
    PreOutPut[3] = X20;
    PreOutPut[4] = X10;
    PreOutPut[5] = kappaStarNew;
    return PreOutPut;
}

//下面是读取数据的算法

double* ReadMicroAndMinus(string MacroType, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Columns, long DisColumns, long DisRows, long EspBoost)
{
    //FIXME最好在微透镜天区里面转，任何一个边都不要超出边界。
    long ColumnsNew = Columns; //这个是错位以后的列数。
    long FileColumn1 = ColumnsNew / EspBoost + DisColumns; //得到初始文件中的列数。
    long FileRow1 = DisRows;
    long RemainNum = ColumnsNew % EspBoost; //得到两列之间的余数。
    double Coeffi = (double) RemainNum / EspBoost;
    vector<double> FileX2SetSub;
    for(long X2SubIndex = 0; X2SubIndex < FileLenRows; X2SubIndex ++ )
    {
        FileX2SetSub.push_back(FileX2Set[X2SubIndex]);
    }
    char FileKappaStar[100];
    sprintf(FileKappaStar,"/disk1/home/shanxk/work/Paper4_MicroLens_Modify/stellar_and_remnant_data_Phi/%3.2f.bin", kappa_star);
    // cout << FileKappaStar << endl;
    /*读取Psi的数据*/
    /*下面是读取数据所在列的下界*/
    double* PsiPerColumnSubFun1 = new double [FileLenRows];
    vector<double> PsiPerColumnSubFunVec1;
    /*数据类型长度*/
    int SizeLen = sizeof(double);
    ifstream PsiFile1;
    if(MacroType == "Minimum")
    {
        PsiFile1.open(FileKappaStar, std::ifstream::binary);
        // PsiFile1.open("/data/pubdata/Data_cxc_to_sxk/sxk/MinimumMicroAndMinus.bin", std::ifstream::binary);
    }
    else if(MacroType == "Maximum")
    {
        PsiFile1.open(FileKappaStar, std::ifstream::binary); 
        // PsiFile1.open("/data/pubdata/Data_cxc_to_sxk/sxk/MaximumMicroAndMinus.bin", std::ifstream::binary); 
    }
    else
    {
        cout << "Unknown Macro Type" << endl;
    }
    if(FileColumn1 <= MicroFileLenColumns -1)
    {

    
        PsiFile1.seekg(SizeLen * (MicroFileLenRows * FileColumn1 + FileRow1), ios::beg); //从第Rows行，第Startpointx这个点开始读数据，其中6是这个文件的数据类型数，第一二列为像的坐标（y,x），第三四列为对应源的坐标(y,x)，第五列为放大率，第六列为引力势。
        if(FileRow1 + FileLenRows <= MicroFileLenRows)
        {
            PsiFile1.read(reinterpret_cast<char *>(PsiPerColumnSubFun1), SizeLen*FileLenRows); //读一整行
            PsiFile1.close ();
        }
        else
        {
            cout << "RowIndex Exceed Micro File" << endl;
            PsiFile1.read(reinterpret_cast<char *>(PsiPerColumnSubFun1), SizeLen*(MicroFileLenRows - FileRow1)); //读一整行
            PsiFile1.close ();
            for(long tmpi = MicroFileLenRows - FileRow1; tmpi < FileLenRows; tmpi ++ )
            {
                PsiPerColumnSubFun1[tmpi] = 0;
            }

        }
        for(long i = 0; i < FileLenRows; i ++ )
        {
            PsiPerColumnSubFunVec1.push_back(PsiPerColumnSubFun1[i]);
        }
        
        tk::spline InterFunc1(FileX2SetSub,PsiPerColumnSubFunVec1);
        
        double* PsiPerColumnSubFunInter1 = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {
            PsiPerColumnSubFunInter1[i] = InterFunc1(X2Set[i]);
        }

        
        


        /*读取Psi的数据*/
        /*下面是读取数据所在行的上界*/
        double* PsiPerColumnSubFun2 = new double [FileLenRows];
        vector<double> PsiPerColumnSubFunVec2;
        /*数据类型长度*/

        ifstream PsiFile2;
        if(MacroType == "Minimum")
        {
            PsiFile2.open(FileKappaStar, std::ifstream::binary);
            // PsiFile2.open("/data/pubdata/Data_cxc_to_sxk/sxk/MinimumMicroAndMinus.bin", std::ifstream::binary);
        }
        else if(MacroType == "Maximum")
        {
            PsiFile2.open(FileKappaStar, std::ifstream::binary);
            // PsiFile2.open("/data/pubdata/Data_cxc_to_sxk/sxk/MaximumMicroAndMinus.bin", std::ifstream::binary); 
        }
        else
        {
            cout << "Unknow Macro Type" << endl;
        }
        if(FileColumn1 < MicroFileLenColumns - 1)
        {
        PsiFile2.seekg(SizeLen * (MicroFileLenRows * ( FileColumn1 + 1 ) + FileRow1 ),ios::beg); //从第Rows行开始读数据 
        }
        else
        {
            PsiFile2.seekg(SizeLen * (MicroFileLenRows * ( FileColumn1 + 0 ) + FileRow1 ),ios::beg); //从第Rows行开始读数据
        }

        if(FileRow1 + FileLenRows <= MicroFileLenRows)
        {
            PsiFile2.read(reinterpret_cast<char *>(PsiPerColumnSubFun2), SizeLen*FileLenRows); //读一整行
            PsiFile2.close ();
        }
        else
        {
            cout << "RowIndex Exceed Micro File" << endl;
            PsiFile2.read(reinterpret_cast<char *>(PsiPerColumnSubFun2), SizeLen*(MicroFileLenRows - FileRow1)); //读一整行
            PsiFile2.close ();
            for(long tmpi = MicroFileLenRows - FileRow1; tmpi < FileLenRows; tmpi ++ )
            {
                PsiPerColumnSubFun2[tmpi] = 0;
            }

        }


        for(long i = 0; i < FileLenRows; i ++ )
        {
            PsiPerColumnSubFunVec2.push_back(PsiPerColumnSubFun2[i]);
        }
        tk::spline InterFunc2(FileX2SetSub,PsiPerColumnSubFunVec2);
        
        // cout << IndexInter << " " << X1SetInter[IndexInter-1] << endl;
        double* PsiPerColumnSubFunInter2 = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {
            PsiPerColumnSubFunInter2[i] = InterFunc2(X2Set[i]);
        }


        /*输出最后插值以后的结果*/
        double* PsiPerColumnSubFunInter = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {

            PsiPerColumnSubFunInter[i] = PsiPerColumnSubFunInter1[i] + (PsiPerColumnSubFunInter2[i] - PsiPerColumnSubFunInter1[i]) * Coeffi;
        }


    


        delete [] PsiPerColumnSubFun1;
        PsiPerColumnSubFun1 = NULL;
        delete [] PsiPerColumnSubFun2;
        PsiPerColumnSubFun2 = NULL;
        delete [] PsiPerColumnSubFunInter1;
        PsiPerColumnSubFunInter1 = NULL; 
        delete [] PsiPerColumnSubFunInter2; 
        PsiPerColumnSubFunInter2 = NULL;
        return PsiPerColumnSubFunInter;
    }

    else //超出范围的部分归零。
    {
        cout << "ColumnIndex Exceed Micro File" << endl;
        double* PsiPerColumnSubFunInter = new double [BoostLenRows]();
        return PsiPerColumnSubFunInter;
    }

    
}



//double* (double kappa, double gamma, double SkyLimit, int LengthX, double SourceX, double SourceY)
double* ReadCreatPsi(string MacroType, double kappa, double gamma, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Column, long DisColumns, long DisRows, double SourceY1, double SourceY2, long EspBoost)
{
    /*计算总的时间延迟，注意这里是无量纲的，返回值是第Column列的所有时间延迟组成的数组。*/
    double* ReadFileArray = ReadMicroAndMinus(MacroType, kappa_star, BoostLenRows, MicroFileLenRows, MicroFileLenColumns, FileLenRows, FileX2Set, X2Set, Column, DisColumns, DisRows, EspBoost);
    double* TotalTimeDelay = new double [BoostLenRows]; 
    for(long RowIndex = 0; RowIndex < BoostLenRows; RowIndex ++ )
    {
        double GeometryDelay = 0.5*(pow((X2Set[Column] - SourceY1),2.) + pow((X2Set[RowIndex] - SourceY2),2.)); //这个是时间延迟的第一项，一般叫做几何项。
        double SmoothDelay = 0.5 * kappa *(X2Set[Column] * X2Set[Column] + X2Set[RowIndex] * X2Set[RowIndex]) 
            + 0.5 * gamma * (X2Set[RowIndex] * X2Set[RowIndex] - X2Set[Column] * X2Set[Column]);
        TotalTimeDelay[RowIndex] = GeometryDelay - SmoothDelay -  ReadFileArray[RowIndex]; //计算总体时间延迟。
    }
    delete [] ReadFileArray;
    ReadFileArray = NULL; 
    
    return TotalTimeDelay;
}




//This is for multi micro lenses, one need to set the \kappa_*.
double* SaddlePreparation4CreatPhiKappaStar(string MacroType, double kappa, double gamma, double kappaStar, double MassPerSolarMass, double LensRedshift, int PrecisionFactor)
{
    // double* PreOutPut = new double[6]; //这个是用来存后面得到的最大值和最小值的。
    double MSun = 2 * pow(10,30);
    double GravityConstant = 6.67 * pow(10,-11);
    double LightSpeed = 3 * pow(10,8);
    double MicroLensCoeffi = 4 * GravityConstant * MassPerSolarMass * MSun * (1 + LensRedshift) / pow(LightSpeed, 3);
    // cout << pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2) << endl;
    double RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)});
    while(RNeed2Cal>900) //这里写的是900，为了留出来200的空隙
    {
        PrecisionFactor -= 1;
        RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}); 
    }
    // cout << "SNR = 60 exceed sky map" << endl;
    cout << "NOW SNR = " << PrecisionFactor << endl;
    double TimeNeed2CalMax = MicroLensCoeffi/2*min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)})*pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2);
    cout << "TimeNeed2CalMax = " << TimeNeed2CalMax << endl;
    double MagnificationR = 1 - kappa + gamma;
    double MagnificationT = 1 - kappa - gamma;

    
    if(MacroType == "Saddle")
    {
        double* PreOutPut = new double[9]; //这个是用来存后面得到的最大值和最小值的。
        double SkyLimitX;
        double SkyLimitY;
        double X20, X10, X2010, X20Prime, X1020, X10Prime, TestX20X10, TestX20X10Max, X10New, X20New;
        double EllipticalA = sqrt(1/MagnificationR/MicroLensCoeffi);
        double EllipticalB = sqrt(-1/MagnificationT/MicroLensCoeffi); //FIXME如果不是magnificationT<0的话，这里需要改。
        double Epsilon1 =35;//50;//30;
        double Epsilon2 = EllipticalB/EllipticalA*Epsilon1; 
        cout << "EllipticalA = " << EllipticalA << " " << "EllipticalB = " << EllipticalB << endl;
        cout << "Epsilon for x20 = " << Epsilon1 << " " << "Epsilon for x10 = " << Epsilon2  << endl;

        X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
        X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
        //X20Prime = EllipticalB/EllipticalA*X2010;
        X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
        X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
        //X10Prime = EllipticalA/EllipticalB*X1020;
        SkyLimitY = X1020;//max({X20Prime, X1020});
        SkyLimitX = X2010;//max({X2010,X10Prime});
        TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
        TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 

        
        cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
        cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
        cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
        

        while((TestX20X10 < 10)&&(Epsilon1>=1)) //防止太小
        {

            Epsilon1 -= 0.5;
            Epsilon2 = EllipticalB/EllipticalA*Epsilon1;
            X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
            X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
            //X20Prime = EllipticalB/EllipticalA*X2010;
            X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
            X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
            //X10Prime = EllipticalA/EllipticalB*X1020;
            SkyLimitY = X1020;//max({X20Prime, X1020});
            SkyLimitX = X2010;//max({X2010,X10Prime});
            TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
            TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 
        }
        cout << "After judge TestX20X10 < 10 and Epsilon1 >= 1" << endl;
        cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
        cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
        cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;

        while(TestX20X10Max > 1000)
        {

            Epsilon1 += 0.5;
            Epsilon2 = EllipticalB/EllipticalA*Epsilon1;
            X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
            X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
            //X20Prime = EllipticalB/EllipticalA*X2010;
            X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
            X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
            //X10Prime = EllipticalA/EllipticalB*X1020;
            SkyLimitY = X1020;//max({X20Prime, X1020});
            SkyLimitX = X2010;//max({X2010,X10Prime});
            TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
            TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 
        }

        cout << "After judge TestX20X10Max > 1000" << endl;
        cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
        cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
        cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
    X10New = X10;//pow(-2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*SkyLimitY*SkyLimitY,0.5);
    X20New = X20;//pow(-2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*SkyLimitX*SkyLimitX,0.5); 

    cout << "Final Epsilon" << endl;
    cout << "Epsilon for X20 = " << Epsilon1 << " Epsilon for X10 = " << Epsilon2 << endl;

    double AreaS = 4 * SkyLimitX * SkyLimitY; //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
    int NStar = kappaStar * AreaS / PI;
    double kappaStarNew = NStar * PI / AreaS;
    cout << "Nstar = " << NStar << endl;
    PreOutPut[0] = X10New;
    PreOutPut[6] = X20New;
    PreOutPut[1] = NStar;
    PreOutPut[2] = TimeNeed2CalMax;
    PreOutPut[3] = SkyLimitX;
    PreOutPut[4] = SkyLimitY;
    PreOutPut[5] = kappaStarNew;
    PreOutPut[7] = Epsilon1;
    PreOutPut[8] = Epsilon2;
    return PreOutPut;
        

        
    }
    else
    {
        cout << "Unknown Macrolensing Type" << endl;
        return 0;
    }
    
    
}



//下面是读取数据的算法

double* SaddleReadMicroAndMinus(double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Columns, long DisColumns, long DisRows, long EspBoost)
{//FileLenRows 是划定的天区行数。FileX1Set是和微透镜+负质量片图分辨率相同的坐标。
    //FIXME最好在微透镜天区里面转，任何一个边都不要超出边界。
    long ColumnsNew = Columns; //这个是错位以后的列数。
    long FileColumn1 = ColumnsNew / EspBoost + DisColumns; //得到初始文件中的列数。
    long RemainNum = ColumnsNew % EspBoost; //得到两列之间的余数。
    double Coeffi = (double) RemainNum / EspBoost;
    vector<double> FileX2SetSub;
    for(long X2SubIndex = 0; X2SubIndex < FileLenRows; X2SubIndex ++ )
    {
        FileX2SetSub.push_back(FileX2Set[X2SubIndex]);
    }
    char FileKappaStar[100];
    sprintf(FileKappaStar,"/disk1/home/shanxk/work/Paper4_MicroLens_Modify/stellar_and_remnant_data_Phi/%3.2f.bin", kappa_star);
    /*读取Psi的数据*/
    /*下面是读取数据所在列的下界*/
    double* PsiPerColumnSubFun1 = new double [FileLenRows];
    vector<double> PsiPerColumnSubFunVec1;
    /*数据类型长度*/
    int SizeLen = sizeof(double);
    ifstream PsiFile1;
    PsiFile1.open(FileKappaStar, std::ifstream::binary);
    if(FileColumn1 <= MicroFileLenColumns -1)
    {

    
        PsiFile1.seekg(SizeLen * (MicroFileLenRows * FileColumn1 + DisRows), ios::beg); //从第Rows行，第Startpointx这个点开始读数据，其中6是这个文件的数据类型数，第一二列为像的坐标（y,x），第三四列为对应源的坐标(y,x)，第五列为放大率，第六列为引力势。
        if(DisRows + FileLenRows <= MicroFileLenRows)
        {
            PsiFile1.read(reinterpret_cast<char *>(PsiPerColumnSubFun1), SizeLen*FileLenRows); //读一整行
            PsiFile1.close ();
        }
        else
        {
            cout << "RowIndex Exceed Micro File" << endl;
            PsiFile1.read(reinterpret_cast<char *>(PsiPerColumnSubFun1), SizeLen*(MicroFileLenRows - DisRows)); //读一整行
            PsiFile1.close ();
            for(long tmpi = MicroFileLenRows - DisRows; tmpi < FileLenRows; tmpi ++ )
            {
                PsiPerColumnSubFun1[tmpi] = 0;
            }

        }
        for(long i = 0; i < FileLenRows; i ++ )
        {
            PsiPerColumnSubFunVec1.push_back(PsiPerColumnSubFun1[i]);
        }
        
        tk::spline InterFunc1(FileX2SetSub,PsiPerColumnSubFunVec1);
        
        double* PsiPerColumnSubFunInter1 = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {
            PsiPerColumnSubFunInter1[i] = InterFunc1(X2Set[i]);
        }

        
        


        /*读取Psi的数据*/
        /*下面是读取数据所在行的上界*/
        double* PsiPerColumnSubFun2 = new double [FileLenRows];
        vector<double> PsiPerColumnSubFunVec2;
        /*数据类型长度*/

        ifstream PsiFile2;
        PsiFile2.open(FileKappaStar, std::ifstream::binary);
        if(FileColumn1 < MicroFileLenColumns - 1)
        {
        PsiFile2.seekg(SizeLen * (MicroFileLenRows * ( FileColumn1 + 1 ) + DisRows ),ios::beg); //从第Rows行开始读数据 
        }
        else
        {
            PsiFile2.seekg(SizeLen * (MicroFileLenRows * ( FileColumn1 + 0 ) + DisRows ),ios::beg); //从第Rows行开始读数据
        }

        if(DisRows + FileLenRows <= MicroFileLenRows)
        {
            PsiFile2.read(reinterpret_cast<char *>(PsiPerColumnSubFun2), SizeLen*FileLenRows); //读一整行
            PsiFile2.close ();
        }
        else
        {
            cout << "RowIndex Exceed Micro File" << endl;
            PsiFile2.read(reinterpret_cast<char *>(PsiPerColumnSubFun2), SizeLen*(MicroFileLenRows - DisRows)); //读一整行
            PsiFile2.close ();
            for(long tmpi = MicroFileLenRows - DisRows; tmpi < FileLenRows; tmpi ++ )
            {
                PsiPerColumnSubFun2[tmpi] = 0;
            }

        }


        for(long i = 0; i < FileLenRows; i ++ )
        {
            PsiPerColumnSubFunVec2.push_back(PsiPerColumnSubFun2[i]);
        }
        tk::spline InterFunc2(FileX2SetSub,PsiPerColumnSubFunVec2);
        
        // cout << IndexInter << " " << X1SetInter[IndexInter-1] << endl;
        double* PsiPerColumnSubFunInter2 = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {
            PsiPerColumnSubFunInter2[i] = InterFunc2(X2Set[i]);
        }


        /*输出最后插值以后的结果*/
        double* PsiPerColumnSubFunInter = new double [BoostLenRows];
        for(long i = 0; i < BoostLenRows; i ++ )
        {

            PsiPerColumnSubFunInter[i] = PsiPerColumnSubFunInter1[i] + (PsiPerColumnSubFunInter2[i] - PsiPerColumnSubFunInter1[i]) * Coeffi;
        }


    


        delete [] PsiPerColumnSubFun1;
        PsiPerColumnSubFun1 = NULL; 
        delete [] PsiPerColumnSubFun2;
        PsiPerColumnSubFun2 = NULL;
        delete [] PsiPerColumnSubFunInter1;
        PsiPerColumnSubFunInter1 = NULL;
        delete [] PsiPerColumnSubFunInter2; 
        PsiPerColumnSubFunInter2 = NULL;
        return PsiPerColumnSubFunInter;
    }

    else //超出范围的部分归零。
    {
        cout << "ColumnIndex Exceed Micro File" << endl;
        double* PsiPerColumnSubFunInter = new double [BoostLenRows]();
        return PsiPerColumnSubFunInter;
    }

    
}



//double* (double kappa, double gamma, double SkyLimit, int LengthX, double SourceX, double SourceY)
double* SaddleReadCreatPsi(double kappa, double gamma, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X1Set, double * X2Set, long Column, long DisColumns, long DisRows, double SourceY1, double SourceY2, long EspBoost)
{
    /*计算总的时间延迟，注意这里是无量纲的，返回值是第Column列的所有时间延迟组成的数组。*/
    double* ReadFileArray = SaddleReadMicroAndMinus(kappa_star, BoostLenRows, MicroFileLenRows, MicroFileLenColumns, FileLenRows, FileX2Set, X2Set, Column, DisColumns, DisRows, EspBoost);
    double* TotalTimeDelay = new double [BoostLenRows]; 
    for(long RowIndex = 0; RowIndex < BoostLenRows; RowIndex ++ )
    {
        double GeometryDelay = 0.5*(pow((X1Set[Column] - SourceY1),2.) + pow((X2Set[RowIndex] - SourceY2),2.)); //这个是时间延迟的第一项，一般叫做几何项。
        double SmoothDelay = 0.5 * kappa *(X1Set[Column] * X1Set[Column] + X2Set[RowIndex] * X2Set[RowIndex]) 
            + 0.5 * gamma * (X2Set[RowIndex] * X2Set[RowIndex] - X1Set[Column] * X1Set[Column]);
        TotalTimeDelay[RowIndex] = GeometryDelay - SmoothDelay -  ReadFileArray[RowIndex]; //计算总体时间延迟。
    }
    delete [] ReadFileArray;
    ReadFileArray = NULL;

    return TotalTimeDelay;
}