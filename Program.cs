using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.Threading.Tasks;
using System.Numerics;


/* Hi! Welcome to my first GitHub upload. I'm a chemical engineering student and
 * I made this program to make some calculations in the lab easier, this program
 * use the gradient descent method to obtain the parameters that better adjust 
 * to the ecuation y=d-(a-d)/(1+c*x^b), using the Least Squares aproach.
 * Because this app is intended to use in a lab the variables for the data are called "Med" (from mediciones)
 */

namespace FinalSolution
{
    class Program
    {

        /*Programa ajustador de parámetros para curva sigmoidea.
         * 
         * Autor: Federico Benelli
         * Versión: 0.0.1
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
            for(int i = 0; i < 4; i++) { Grad[i] = 0.0; }
            

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
            double gamma=0;

            
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
                est[i] = d - (a - d) / (1 + c*exp);
            }
            return est;
        }

        //Calcula el error cuadrado total
        public static double CuadErr(double[] y, double[] est)
        {
            double errCuad = 0;

            for (int i=0;i<y.Length;i++)
            {
                errCuad += Math.Pow(y[i] - est[i], 2);
            }

            return errCuad;
        }

        static void Main(string[] args)
        {
            
            
            Console.WriteLine("Bienvenido a Calculador de Amilosa v0.01");
            Console.WriteLine("A continuación deberá agregar los valores que se le soliciten, cuando halla agregado el valor que considera suficiente presione enter sin agregar nada.")
            bool h=True;
            double[] vectPar={-10,0,32,2};
            double[] vectVal={};
            double[] x={};
            int addVal=0;
            
            //Here the user adds the values of x and y
            while(h)
            {
                double[] temp={0};
                if(temp[0]!=null)
                {
                    Console.Write("Volumen {0}: ",addVal);
                    temp[0]=Console.Read();
                    x=x.Concat(temp).ToArray();
                    Console.WriteLine();
                    Console.Write("Medición {0}: ",addVal);
                    temp[0]=Console.Read();
                    vectVal=vectVal.Concat(temp).ToArray();
                }
                else{h=False;}
                addVal+=1;
            }
            
            double[] prueba = { 0 };
            double[] est;
            double error;
            double[] estIni;
            double errorIni;
            double[] vectParini = new double[vectPar.Length];
            for (int i = 0; i < vectPar.Length; i++) { vectParini[i] = vectPar[i]+2; }
            int h = 0;
            double gamma = 0.01;


            //Iterations of the gradient descent
            while (h < 5000)
            {
                
                est = Sig(x, vectPar);
                double[] vectGrad = Grad(x, vectPar, vectVal);
                double[] vectGradini = Grad(x, vectParini, vectVal);
                double[] rest =new double[vectGrad.Length];
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
                
                est = Sig(x, vectPar);
                estIni = Sig(x, vectParini);
                error = CuadErr(vectVal, est);
                errorIni = CuadErr(vectVal, estIni);

                for (int i = 0; i < vectPar.Length; i++) { vectParini[i] = vectPar[i]; }
                
                h += 1;
            }

            Console.WriteLine("Los parámetros son: ");
            foreach(double i in vectPar) { Console.Write("{0};  ", i); }
            Console.Read();
        }
    }
}
