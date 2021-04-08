using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using static System.Console;
using System.Diagnostics;

namespace Auros1stProject2_2_MyAlgorithm
{
    class Program
    {
        static void Main(string[] args)
        {
            //
            // "SiO2 1000nm_on_Si.dat" 파일 로딩 후
            // 측정 스펙트럼 데이터를 alpha, beta 로 변환한다.
            //
            // 2021.03.24 이지원.
            //
            #region psi, delta -> alpha, beta

            List<string> MeasurementSpectrumData = new List<string>();  // 측정 스펙트럼 데이터 저장할 배열. (한 줄씩 저장)
            string[] SingleLineData;                                    // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            // "SiO2 2nm_on_Si.dat" 파일 읽기. (한 줄씩)
            MeasurementSpectrumData.AddRange(File.ReadAllLines("SiO2 1000nm_on_Si.dat"));

            // 무의미한 공백 행을 제거한다.
            int lenSpectrumData = MeasurementSpectrumData.Count;
            string Blank = "";
            for (int i = 0; i < lenSpectrumData; i++)
            {
                if (MeasurementSpectrumData[0] == Blank)
                    MeasurementSpectrumData.RemoveAt(0);
                else
                    break;
            }

            // wavelength : 350 ~ 980(nm)인 측정 스펙트럼 데이터를 담을 리스트 선언.
            List<double> wavelength_exp = new List<double>();   // 파장 데이터 리스트.
            List<double> AOI_exp = new List<double>();          // 입사각 데이터 리스트.
            List<double> alpha_exp = new List<double>();        // Psi 데이터 리스트.
            List<double> beta_exp = new List<double>();         // Delta 데이터 리스트.

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하기 위해 반복문을 1부터 시작한다.
            int StartIndex = 1;
            int LenData = MeasurementSpectrumData.Count;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = MeasurementSpectrumData[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 파장이 350 ~ 980(nm) 이내인 데이터만 저장한다.
                if (Convert.ToDouble(SingleLineData[0]) >= 350.0 &&
                    Convert.ToDouble(SingleLineData[0]) <= 980.0)
                {
                    // 각 컬럼에 해당하는 데이터를 저장한다.
                    wavelength_exp.Add(Double.Parse(SingleLineData[0]));
                    AOI_exp.Add(Double.Parse(SingleLineData[1]));
                    alpha_exp.Add(Double.Parse(SingleLineData[2]));
                    beta_exp.Add(Double.Parse(SingleLineData[3]));
                }
            }

            // psi, delta -> alpha, beta 변환.

            // degree, radian 변환 인라인 함수 정의.
            double degree2radian(double angle) => ((angle * (PI)) / 180.0);
            //double radian2degree(double angle) => (angle * (180.0 / PI));

            // Polarizer offset 각도. (45도)
            double PolarizerRadian = degree2radian(45.0);

            // psi, delta 데이터를 alpha, beta 로 변환한다.
            LenData = wavelength_exp.Count;
            for (int i = 0; i < LenData; i++)
            {
                // psi, delta 값을 radian 으로 변환한다.
                double PsiRadian = degree2radian(alpha_exp[i]);
                double DeltaRadian = degree2radian(beta_exp[i]);

                // psi, delta 데이터를 alpha, beta 로 갱신한다.
                alpha_exp[i] = (
                    (Pow(Tan(PsiRadian), 2.0) - Pow(Tan(PolarizerRadian), 2.0))
                    / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
                beta_exp[i] = (
                    (2.0 * Tan(PsiRadian) * Tan(PolarizerRadian) * Cos(DeltaRadian))
                    / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
            }

            #endregion

            //
            // "Si_new.txt", "SiO2_new.txt" 파일 물성값 로딩.
            //
            // 2021.03.24 이지원.
            //
            #region MyRegion

            // "Si_new.txt" 파일 읽기.
            string[] Si_new = File.ReadAllLines("Si_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = Si_new.Length - 1;
            double[] wavelength_Si = new double[LenData];
            double[] n_Si = new double[LenData];
            double[] k_Si = new double[LenData];

            // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = Si_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = Si_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_Si[i - 1] = Double.Parse(SingleLineData[0]);
                n_Si[i - 1] = Double.Parse(SingleLineData[1]);
                k_Si[i - 1] = Double.Parse(SingleLineData[2]);
            }


            // "SiO2_new.txt" 파일 읽기.
            string[] SiO2_new = File.ReadAllLines("SiO2_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = SiO2_new.Length - 1;
            double[] wavelength_SiO2 = new double[LenData];
            double[] n_SiO2 = new double[LenData];
            double[] k_SiO2 = new double[LenData];

            // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = SiO2_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiO2_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiO2[i - 1] = Double.Parse(SingleLineData[0]);
                n_SiO2[i - 1] = Double.Parse(SingleLineData[1]);
                k_SiO2[i - 1] = Double.Parse(SingleLineData[2]);
            }

            #region Si_new, SiO2_new 데이터 출력 (Test)
            /*LenData = wavelength_Si.Length;
            for (int i = 0; i < LenData; i++)
                WriteLine(wavelength_Si[i] + "\t" + n_Si[i] + "\t" + k_Si[i]);
            WriteLine("============================================");
            for (int i = 0; i < LenData; i++)
                WriteLine(wavelength_SiO2[i] + "\t" + n_SiO2[i] + "\t" + k_SiO2[i]);*/
            #endregion
            #endregion

            //
            // "Si_new.txt", "SiO2_new.txt" 의 n, k 를 사용하여
            // 각 계면에서의 반사, 투과계수를 계산한다.
            //
            // 2021.03.24 이지원.
            //
            #region 각 계면에서의 반사, 투과계수 계산

            LenData = wavelength_Si.Length;

            // 반사계수를 담을 배열.
            Complex[] r12p = new Complex[LenData],
                      r12s = new Complex[LenData],
                      r01p = new Complex[LenData],
                      r01s = new Complex[LenData];
            // 투과계수를 담을 배열.
            Complex[] t12p = new Complex[LenData],
                      t12s = new Complex[LenData],
                      t01p = new Complex[LenData],
                      t01s = new Complex[LenData];

            double AOI_air = degree2radian(65.0);   // 입사각. (라디안) 
            Complex N_air = new Complex(1.0, 0);    // 공기의 굴절률.

            // 반사, 투과계수를 계산한다.
            for (int i = 0; i < LenData; i++)
            {
                // 파장에 대한 물질의 복소굴절률을 구한다.
                Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);
                Complex N_Si = new Complex(n_Si[i], -k_Si[i]);

                // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Complex Sintheta_j = new Complex(Sin(AOI_air), 0);
                Complex Costheta_j = new Complex(Cos(AOI_air), 0);
                Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                Complex theta_k = Complex.Asin(Sintheta_k);
                // air, SiO2 경계면에서의 굴절각.
                Complex Costheta_k = Complex.Cos(theta_k);

                // air, SiO2 경계면에서의 반사계수를 구한다.
                r01p[i] = ((N_SiO2 * Costheta_j) - (N_air * Costheta_k)) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                r01s[i] = ((N_air * Costheta_j) - (N_SiO2 * Costheta_k)) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                // air, SiO2 경계면에서의 투과계수를 구한다.
                t01p[i] = (N_air * Costheta_j * 2.0) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                t01s[i] = (N_air * Costheta_j * 2.0) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                // SiO2, Si 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Sintheta_j = Complex.Sin(theta_k);
                Costheta_j = Complex.Cos(theta_k);
                Sintheta_k = (N_SiO2 / N_Si) * Sintheta_j;
                theta_k = Complex.Asin(Sintheta_k);         // SiO2, Si 경계면에서의 굴절각.
                Costheta_k = Complex.Cos(theta_k);

                // SiO2, Si 경계면에서의 반사계수를 구한다.
                r12p[i] = ((N_Si * Costheta_j) - (N_SiO2 * Costheta_k)) /
                             ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                r12s[i] = ((N_SiO2 * Costheta_j) - (N_Si * Costheta_k)) /
                             ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));

                // SiO2, Si 경계면에서의 투과계수를 구한다.
                t12p[i] = (N_SiO2 * Costheta_j * 2.0) /
                             ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                t12s[i] = (N_SiO2 * Costheta_j * 2.0) /
                             ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));
            }

            #region 위에서 구한 반사, 투과계수 출력 (Test)
            /*WriteLine("====== air, SiO2 경계 ======");
            for (int i = 0; i < LenData; i++)
            {
                WriteLine(
                    r01p[i] + " " +
                    r01s[i] + " " +
                    t01p[i] + " " +
                    t01s[i]);
            }
            WriteLine("====== SiO2, Si 경계 ======");
            for (int i = 0; i < LenData; i++)
            {
                WriteLine(
                    r12p[i] + " " +
                    r12s[i] + " " +
                    t12p[i] + " " +
                    t12s[i]);
            }*/
            #endregion

            #endregion


            // MSE 와 그에 대응하는 두께를 저장할 리스트 선언.
            double nowMSE = 0.0,
                   preMSE = 0.0;
            double nowThickness = 0.0,
                   preThickness = 0.0;

            double d0 = 800.0;                      // 탐색을 시작할 두께.
            const double ThicknessRecipe = 1000.0;  // 레시피 두께.



            double StepSize = 0.9;  // step size.
            double scale = 1.1;     // step size 스케일 값.
            int IsLeft = 1;         // d0 위치 판단 변수.

            #region 방향을 판단한다.

            // d0 가 ThicknessRecipe 보다 오른쪽에 있으면 탐색 방향을 왼쪽으로 설정한다.
            if (ThicknessRecipe - d0 < 0)
            {
                StepSize = -StepSize;
                IsLeft = 0;
            }
            #endregion
            #region d0 에 대한 첫번째 두께, MSE 를 계산한다.

            #region 총 반사계수를 계산한다.

            Complex[] Rp = new Complex[LenData],
                      Rs = new Complex[LenData];
            for (int i = 0; i < LenData; i++)
            {
                // SiO2의 복소 굴절률.
                Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                Complex theta_k = Complex.Asin(Sintheta_k);         // air, SiO2 경계면에서의 굴절각.
                Complex Costheta_k = Complex.Cos(theta_k);

                // 위상 두께를 구한다.
                Complex PhaseThickness = ((double)d0 * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                //WriteLine(PhaseThickness);

                // 총 반사계수를 구한다.
                Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                Rp[i] = (r01p[i] + r12p[i] * E) /
                        (1 + r01p[i] * r12p[i] * E);

                Rs[i] = (r01s[i] + r12s[i] * E) /
                        (1 + r01s[i] * r12s[i] * E);

            }

            #endregion
            #region 총 반사계수로부터 alpha, beta 를 도출한다.

            // alpha, beta 이론값을 담을 배열 선언.
            double[] alpha_cal = new double[LenData],
                     beta_cal = new double[LenData];

            // Polarizer 오프셋 각.
            double polarizerAngle = degree2radian(45.0);

            for (int i = 0; i < LenData; i++)
            {
                // 총 반사계수비. (복소반사계수비)
                Complex rho = Rp[i] / Rs[i];

                // Psi, Delta.
                double Psi = Atan(rho.Magnitude);
                double Delta = rho.Phase;


                alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                       (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                       (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
            }
            #region alpha, beta 이론값 출력 (Test)

            /*for (int i = 0; i < LenData; i++)
                WriteLine(alpha_cal[i] + " " 
                + beta_cal[i]);*/

            #endregion
            #endregion
            #region 첫번째 MSE 를 계산하고 MSE 와 그 때의 두께를 저장한다.

            double sum = 0;
            for (int i = 0; i < LenData; i++)
            {
                double difference_MSE =
                     Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                     Pow((beta_exp[i] - beta_cal[i]), 2.0);
                sum += difference_MSE;

            }

            preMSE = sum / LenData;   // 두께에 대해 구해진 MSE 값 저장.
            preThickness = d0;        // 현재 MSE 에서의 두께 저장.

            #endregion

            #endregion

            // 현재 두께를 갱신한다.
            nowThickness = d0 + StepSize;

            int times = 0;
            while (true)
            {
                ++times;
                #region 총 반사계수를 계산한다.

                // 총 반사계수를 저장할 배열 선언.
                Rp = new Complex[LenData];
                Rs = new Complex[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // SiO2의 복소 굴절률.
                    Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                    // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                    Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                    Complex theta_k = Complex.Asin(Sintheta_k);         // air, SiO2 경계면에서의 굴절각.
                    Complex Costheta_k = Complex.Cos(theta_k);

                    // 위상 두께를 구한다.
                    Complex PhaseThickness = ((double)nowThickness * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                    //WriteLine(PhaseThickness);

                    // 총 반사계수를 구한다.
                    Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                    Rp[i] = (r01p[i] + r12p[i] * E) /
                            (1 + r01p[i] * r12p[i] * E);

                    Rs[i] = (r01s[i] + r12s[i] * E) /
                            (1 + r01s[i] * r12s[i] * E);

                }
                #endregion
                #region 총 반사계수로부터 alpha, beta 도출.

                // alpha, beta 이론값을 담을 배열 선언.
                alpha_cal = new double[LenData];
                beta_cal = new double[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // 총 반사계수비. (복소반사계수비)
                    Complex rho = Rp[i] / Rs[i];

                    // Psi, Delta.
                    double Psi = Atan(rho.Magnitude);
                    double Delta = rho.Phase;


                    alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                    beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                }

                #region alpha, beta 이론값 출력 (Test)

                /*for (int i = 0; i < LenData; i++)
                    WriteLine(alpha_cal[i] + " " 
                    + beta_cal[i]);*/

                #endregion
                #endregion
                #region MSE 를 계산한다.

                sum = 0;
                for (int i = 0; i < LenData; i++)
                {
                    double difference_MSE =
                         Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                         Pow((beta_exp[i] - beta_cal[i]), 2.0);
                    sum += difference_MSE;

                }

                // 현재 MSE 값을 구한다.
                nowMSE = sum / LenData;
                #endregion
                #region 이전 MSE 와의 기울기를 계산하여 기울기 부호 변화 시 while 문을 탈출한다.

                // 두께 변화량, MSE 변화량, 기울기를 계산한다.
                double dy = nowMSE - preMSE;
                double dx = nowThickness - preThickness;
                double gradient = dy / dx;

                // WriteLine($"{times} \t {nowMSE} \t {nowThickness}");
                //WriteLine(nowMSE);
                WriteLine(nowThickness);

                // 이전 두께를 갱신한다.
                preThickness = nowThickness;

                // 기울기를 판단한다.
                switch (IsLeft)
                {
                    // d0 가 recipe 두께보다 오른쪽에 있을 때.
                    case 0:
                        {
                            // 오르막일 때.
                            if (gradient < 0)
                            {
                                // step size 를 늘려준다.
                                StepSize *= scale;
                                nowThickness += StepSize;
                            }
                            // 내리막일 때.
                            else
                                goto MomentumStartingPoint;
                        }
                        break;

                    // d0 가 recipe 두께보다 왼쪽에 있을 때.
                    case 1:
                        {
                            // 오르막일 때.
                            if (gradient > 0)
                            {
                                // step size 를 늘려준다.
                                StepSize *= scale;
                                nowThickness += StepSize;
                            }
                            else
                                goto MomentumStartingPoint;
                        }
                        break;

                    default:
                        break;
                }

                preMSE = nowMSE;       // 두께에 대해 구해진 MSE 값을 저장한다.

                #endregion
            }

            // 위의 while 문에서 내리막을 찾았을 때 이 지점으로 오게 된다.
            MomentumStartingPoint:;

            // nowThickness를 갱신한다.
            nowThickness += StepSize;
            StepSize = 0.9;

            // 모멘텀을 위한 변수를 선언한다.
            double mt = 0.0, vt = 0.0, mt_biascorr = 0.0, vt_biascorr = 0.0;
            double Epsilon = Pow(10, -8);

            int cnt = 1;    // 회차 변수.

            // global minimum 과 그 때의 두께 값 변수.
            double globalMinimum, d_sol;

            // 내리막을 찾은 후 Adam 알고리즘을 적용한다.
            while (true)
            {
                ++times;
                #region 총 반사계수를 계산한다.

                // 총 반사계수를 저장할 배열 선언.
                Rp = new Complex[LenData];
                Rs = new Complex[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // SiO2의 복소 굴절률.
                    Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                    // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                    Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                    Complex theta_k = Complex.Asin(Sintheta_k);         // air, SiO2 경계면에서의 굴절각.
                    Complex Costheta_k = Complex.Cos(theta_k);

                    // 위상 두께를 구한다.
                    Complex PhaseThickness = ((double)nowThickness * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                    //WriteLine(PhaseThickness);

                    // 총 반사계수를 구한다.
                    Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                    Rp[i] = (r01p[i] + r12p[i] * E) /
                            (1 + r01p[i] * r12p[i] * E);

                    Rs[i] = (r01s[i] + r12s[i] * E) /
                            (1 + r01s[i] * r12s[i] * E);

                }
                #endregion
                #region 총 반사계수로부터 alpha, beta 도출.

                // alpha, beta 이론값을 담을 배열 선언.
                alpha_cal = new double[LenData];
                beta_cal = new double[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // 총 반사계수비. (복소반사계수비)
                    Complex rho = Rp[i] / Rs[i];

                    // Psi, Delta.
                    double Psi = Atan(rho.Magnitude);
                    double Delta = rho.Phase;


                    alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                    beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                }

                #region alpha, beta 이론값 출력 (Test)

                /*for (int i = 0; i < LenData; i++)
                    WriteLine(alpha_cal[i] + " " 
                    + beta_cal[i]);*/

                #endregion
                #endregion
                #region MSE 를 계산한다.

                sum = 0;
                for (int i = 0; i < LenData; i++)
                {
                    double difference_MSE =
                         Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                         Pow((beta_exp[i] - beta_cal[i]), 2.0);
                    sum += difference_MSE;

                }

                // 현재 MSE 값을 구한다.
                nowMSE = sum / LenData;
                #endregion
                //WriteLine($"{times} \t {nowMSE} \t {nowThickness}");
                //WriteLine(nowMSE);
                WriteLine(nowThickness);

                #region 최적화 완료 조건을 확인한다.

                // MSE 값의 변화가 0.00001(10^-6) 보다 작으면 최적화 수행을 중단한다.
                if (Abs(nowMSE - preMSE) <= 0.000001)
                {
                    // global minimum은 둘중 낮은 MSE를 선택한다.
                    globalMinimum = (nowMSE < preMSE) ? nowMSE : preMSE;

                    // global minimum에 해당하는 두께를 최종 두께인 d_sol로 선정한다.
                    d_sol = (nowMSE == globalMinimum) ? nowThickness : preThickness;

                    goto FindDSol;
                }
                #endregion
                #region Adam 알고리즘을 적용한다.
                // 두께 변화량, MSE 변화량, 기울기를 계산한다.
                double dy = nowMSE - preMSE;
                double dx = nowThickness - preThickness;
                double gradient = dy / dx;
                // double degreeOfGradient = radian2degree(Atan(gradient)) + 90.0;

                // 이전 두께를 갱신한다.
                preThickness = nowThickness;
                preMSE = nowMSE;

                // 아담 알고리즘을 사용한다.
                mt = 0.9 * mt + (1 - 0.9) * gradient;
                vt = 0.999 * vt + (1 - 0.999) * Pow(gradient, 2);
                mt_biascorr = mt / (1 - Pow(0.9, cnt));
                vt_biascorr = vt / (1 - Pow(0.999, cnt));
                nowThickness = preThickness - (StepSize * mt_biascorr) / (Sqrt(vt_biascorr) + Epsilon);

                // 회차를 증가시킨다.
                cnt++;
                #endregion
            }
            FindDSol: WriteLine($"회차 : {times},    global minimum (MSE) : {globalMinimum},   최종 두께 : {d_sol}nm");

        }
    }
}
