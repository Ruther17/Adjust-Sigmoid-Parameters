using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.Threading.Tasks;

namespace Adjust_Sigmoid_Parameters
{
    class Program
    {

        /*Programa ajustador de parámetros para curva sigmoidea.
         * 
         * Autor: Federico Benelli; DNI: 37.300.407
         *
         */


        //Calcula derivadas y devuelve gradiente
        public static double[] Grad(double[] x, double[] par, double[] V)

        {
            double a = par[0];
            double b = par[1];
            double c = par[2];
            double d = par[3];

            double temp;
            int h = 0;
            double[] Grad = new double[4];
            for (int i = 0; i < 4; i++) { Grad[i] = 0.0; }

            foreach (double i in x)
            {
                temp = Math.Pow(i, b);
                Grad[0] += -2 * (-(a - d) / (c * temp + 1) + d - V[h]) / (c * temp + 1);
                Grad[1] += 2 * c * (a - d) * Math.Log(i, Math.E) * temp * (-(a - d) / (c * temp + 1) + d - V[h]) / Math.Pow((c * temp + 1), 2);
                Grad[2] += 2 * (a - d) * temp * (-(a - d) / (c * temp + 1) + d - V[h]) / Math.Pow((c * temp + 1), 2);
                Grad[3] += 2 * (1 / (c * temp + 1) + 1) * (-(a - d) / (c * temp + 1) + d - V[h]);
                h += 1;
            }

            return Grad;
        }

        //Calculo del módulo
        public static double Mod(double[] vect)
        {
            double mod = 0;
            foreach (double i in vect) { mod += Math.Pow(i, 2); }
            mod = (Math.Sqrt(mod));
            return mod;
        }

        //Calculo del gamma
        public static double Gamma(double[] Pinf, double[] Psup, double[] gradInf, double[] gradSup, double mod)
        {
            double gamma = 0;


            for (int i = 0; i < 4; i++)
            {
                gamma += (Psup[i] - Pinf[i]) * (gradSup[i] - gradInf[i]);
            }
            if (mod != 0)
            { gamma = gamma / Math.Pow(mod, 2); }


            return gamma;
        }

        //Suma de matrices
        public static double[,] SumMat(double[,] izq, double[,] der)
        {
            double[,] sum = new double[izq.GetLength(0), der.GetLength(1)];
            for (int i = 0; i < izq.GetLength(0); i++)
            {
                for (int j = 0; j < der.GetLength(1); j++)
                {
                    sum[i, j] = izq[i, j] + der[i, j];
                }
            }
            return sum;
        }

        //Cálculo del volumen equivalente  ----> Corregir las potencias del volumen equivalente. Es lo que da NaN (la primera es la más complicada)
        public static double VolEq(double[] vol, double[] vectPara)
        {
            //Calculo el máximo de la derivada para hallar el punto de inflexión

            double vol1 = 1;
            double vol2 = 2;

            double x = 0.01;
            double a = vectPara[0];
            double b = vectPara[1];
            double c = vectPara[2];
            double d = vectPara[3];
            while (vol2 > vol1)
            {
                vol1 = b * c * (a - d) * Math.Pow(x, b - 1) / Math.Pow((c * Math.Pow(x, b) + 1), 2);
                x += 0.005;
                vol2 = b * c * (a - d) * Math.Pow(x, b - 1) / Math.Pow((c * Math.Pow(x, b) + 1), 2);
            }

            double[] est = Sig(vol, vectPara);
            double P = (est[0] + est[1] + est[2] + est[3] + est[4]) / 5;
            double volEq;
            double l = x;
            volEq = -(Math.Pow(l, 1 - b)) * (d * Math.Pow((c * Math.Pow(l, b) + 1), 2) - P * Math.Pow((c * Math.Pow(l, b) + 1), 2) - (a - d) * (c * Math.Pow(l, b + 1)) - b * c * (a - d) * Math.Pow(l, b)) / (b * c * (a - d));
            return volEq;
        }

        //Resta de matrices
        public static double[,] RestMat(double[,] izq, double[,] der)
        {
            double[,] res = new double[izq.GetLength(0), der.GetLength(1)];
            for (int i = 0; i < izq.GetLength(0); i++)
            {
                for (int j = 0; j < der.GetLength(1); j++)
                {
                    res[i, j] = izq[i, j] - der[i, j];
                }
            }
            return res;
        }

        //Multiplicar matrices
        public static double[,] MultMat(double[,] matIzq, double[,] matDer)
        {
            double[,] matRes = new double[matIzq.GetLength(0), matDer.GetLength(1)];
            for (int i = 0; i < matRes.GetLength(0); i++)
            {
                for (int j = 0; j < matRes.GetLength(1); j++)
                {
                    int r = 0;
                    while (r < matDer.GetLength(0)) { matRes[i, j] += matIzq[i, r] * matDer[r, i]; r += 1; }
                }
            }
            return matRes;
        }

        //Calcula la sigmoidea en los puntos que corresponden a las mediciones y devuelve la lista de valores estimados    
        public static double[] Sig(double[] x, double[] P)
        {
            double exp;

            double[] est = new double[x.Length];
            double a = P[0];
            double b = P[1];
            double c = P[2];
            double d = P[3];
            for (int i = 0; i < x.Length; i++)
            {
                exp = Math.Pow(x[i], b);
                est[i] = d - (a - d) / (1 + c * exp);
            }
            return est;
        }

        //Calcula el error cuadrado total
        public static double CuadErr(double[] y, double[] est)
        {
            double errCuad = 0;

            for (int i = 0; i < y.Length; i++)
            {
                errCuad += Math.Pow(y[i] - est[i], 2);
            }

            return errCuad;
        }

        //Método de gradiente decreciente
        public static void Min(ref double[] x, ref double[] vectPar, ref double[] vectVal)
        {
            double[] est;
            double error;
            double[] estIni;
            double errorIni;
            double[] vectParini = new double[vectPar.Length];
            for (int i = 0; i < vectPar.Length; i++) { vectParini[i] = vectPar[i] + 2; }
            int c = 0;
            double gamma = 0.001;

            bool h = true;

            while (h)
            {

                est = Sig(x, vectPar);
                double[] vectGrad = Grad(x, vectPar, vectVal);
                double[] vectGradini = Grad(x, vectParini, vectVal);
                double[] rest = new double[vectGrad.Length];


                for (int i = 0; i < vectPar.Length; i++)
                {
                    vectParini[i] = vectPar[i];
                }

                for (int i = 0; i < vectPar.Length; i++)
                {
                    vectPar[i] = vectPar[i] - gamma * vectGrad[i];
                }

                vectGrad = Grad(x, vectPar, vectVal);

                for (int i = 0; i < vectPar.Length; i++)
                {
                    rest[i] = vectGrad[i] - vectGradini[i];
                }

                double mod = Mod(rest);
                if (0.001 / Mod(vectGrad) != Double.NaN) { gamma = 0.001 / Mod(vectGrad); }


                est = Sig(x, vectPar);
                estIni = Sig(x, vectParini);
                error = CuadErr(vectVal, est) / vectVal.Length;
                errorIni = CuadErr(vectVal, estIni) / vectVal.Length;

                //Console.WriteLine(error);                                                                          //Seguimiento iteración por iteración

                for (int i = 0; i < vectPar.Length; i++) { vectParini[i] = vectPar[i]; }

                if (Math.Round(error, 1) == Math.Round(errorIni, 1) & error / vectVal.Length > 1000)
                {
                    gamma = 0.5;
                }
                if (Math.Round(error, 10) == Math.Round(errorIni, 10))
                {
                    gamma = 0.001 / mod;
                }
                if (Math.Round(error, 5) == Math.Round(errorIni, 5) & error / vectVal.Length < 10) { break; }

            }
        }

        static void Main(string[] args)
        {

            Console.WriteLine("Bienvenido a Calculador de Amilosa v0.01");
            Console.WriteLine("A continuación deberá agregar los valores que se le soliciten, cuando halla agregado una cantidad de valores que considere suficiente escriba '000'.");
            Console.Write("Inserte el valor de la medición a Volumen cero: ");

            double[] vectPar = { 10, 1, 8, 50 };
            double[] vectVal = {  /*50, 50, 50, 52, 52, 53, 53, 54, 54, 60, 66, 72, 78, 80, 82, 84, 85, 86, 87, 87*/ };                     //valores comentados son para pruebas
            double[] x = { /* 0.0001, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8 */};

            int addVal = 1;
            var h = true;
            var tempo = "";
            double temp;
            //Pedir datos

            /**/
            while (h)
            {
                Console.Write("Volumen {0}: ", addVal);
                tempo = Console.ReadLine();
                tempo = tempo.Replace(".", ",");
                temp = Convert.ToDouble(tempo);
                double[] val = { temp };
                if (temp == 000) { break; }
                x = x.Concat(val).ToArray();

                Console.Write("Medición {0}: ", addVal);
                tempo = Console.ReadLine();
                tempo = tempo.Replace(".", ",");
                temp = Convert.ToDouble(tempo);
                val[0] = temp;
                if (temp == 000) { break; }
                vectVal = vectVal.Concat(val).ToArray();

                addVal++;
            }

            /**/



            //Iteración para ajustar parámetros 

            Min(ref x, ref vectPar, ref vectVal);

            Console.WriteLine("Los parámetros son: ");
            foreach (double i in vectPar) { Console.Write("{0};  ", i); }
            double [] est = Sig(x, vectPar);
            double error = CuadErr(vectVal, est) / vectVal.Length;
            Console.WriteLine();
            Console.WriteLine("El error cuadrado es de: {0}", error);
            Console.WriteLine();
            Console.WriteLine("El volumen equivalente es: {0}", VolEq(x, vectPar));
            Console.Read();
        }
    }
}
