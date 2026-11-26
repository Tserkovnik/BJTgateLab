using System;
using System.Collections.Generic;

namespace BJTgateLab
{
    class Solver
    {
        public int Iters;
        public double FinalNorm;
        
        const double Is = 5.47e-12, Vt = 0.025, Gmin = 1e-8;

        readonly Netlist nl;
        readonly int N, M;                 // N узлов (с землёй), M неизвестных
        readonly int[] branch;             // доп. неизвестная (ток) для каждого V
        readonly Dictionary<int, double[]> prev = new Dictionary<int, double[]>(); // [vbe,vbc] на Q

        public Solver(Netlist nl)
        {
            this.nl = nl;
            N = nl.NodeNames.Count;
            int nV = 0;
            branch = new int[nl.Devices.Count];
            for (int i = 0; i < nl.Devices.Count; i++)
                if (nl.Devices[i].Kind == Kind.V) branch[i] = (N - 1) + nV++;
            M = (N - 1) + nV;
        }

        static int idx(int node) => node == 0 ? -1 : node - 1;              // земля -> -1
        static double Volt(double[] v, int node) => node == 0 ? 0 : v[node - 1];

        // проводимость между узлами a,b
        static void g2(double[,] J, int a, int b, double g)
        {
            if (a >= 0) J[a, a] += g;
            if (b >= 0) J[b, b] += g;
            if (a >= 0 && b >= 0) { J[a, b] -= g; J[b, a] -= g; }
        }

        // управляемый источник: ток между f,t зависит от напряжения ctrl
        static void cross(double[,] J, int ctrl, int f, int t, double g)
        {
            if (ctrl >= 0) J[ctrl, ctrl] += g;
            if (ctrl >= 0 && t >= 0) J[ctrl, t] -= g;
            if (f >= 0 && ctrl >= 0) J[f, ctrl] -= g;
            if (f >= 0 && t >= 0) J[f, t] += g;
        }

        // обрезка скачка напряжения перехода
        static double pnLim(double v1, double v0, double vcrit)
        {
            if (v0 < 0 && v1 > vcrit) return vcrit;
            if (v0 >= 0 && v0 < vcrit && v1 > vcrit) return vcrit;
            double dv = v1 - v0;
            if (v1 > vcrit && Math.Abs(dv) > 2 * Vt)
            {
                double arg = dv / Vt;
                return arg + 1 < 0 ? vcrit : v0 + Vt * Math.Log(1 + arg);
            }
            return v1;
        }

        static double MaxAbs(double[] v)
        {
            double m = 0;
            foreach (var t in v) m = Math.Max(m, Math.Abs(t));
            return m;
        }

        // метод исключения Гаусса с выбором главного элемента
        static double[] Gauss(double[,] a, double[] b)
        {
            int n = b.Length;
            for (int c = 0; c < n; c++)
            {
                int p = c;
                for (int r = c + 1; r < n; r++)
                    if (Math.Abs(a[r, c]) > Math.Abs(a[p, c])) p = r;
                for (int j = 0; j < n; j++) { var t = a[c, j]; a[c, j] = a[p, j]; a[p, j] = t; }
                { var t = b[c]; b[c] = b[p]; b[p] = t; }
                if (Math.Abs(a[c, c]) < 1e-30) continue;
                for (int r = c + 1; r < n; r++)
                {
                    double m = a[r, c] / a[c, c];
                    for (int j = c; j < n; j++) a[r, j] -= m * a[c, j];
                    b[r] -= m * b[c];
                }
            }
            var x = new double[n];
            for (int r = n - 1; r >= 0; r--)
            {
                double s = b[r];
                for (int j = r + 1; j < n; j++) s -= a[r, j] * x[j];
                x[r] = Math.Abs(a[r, r]) < 1e-30 ? 0 : s / a[r, r];
            }
            return x;
        }

        // истинные токи Q (Эберс-Молл) в точке v
        void TrueBJT(Device d, double[] v, double gadd, out double Ic, out double Ib, out double Ie)
        {
            double aF = d.Value / (d.Value + 1);
            double bR = 0.1 * d.Value;
            double aR = Math.Max(0.5, bR / (bR + 1));
            double vbe = Volt(v, d.B) - Volt(v, d.C), vbc = Volt(v, d.B) - Volt(v, d.A);
            double Ebe = Is * (Math.Exp(vbe / Vt) - 1);
            double Ebc = Is * (Math.Exp(vbc / Vt) - 1);
            Ic = aF * Ebe - Ebc;
            Ib = (1 - aF) * Ebe + (1 - aR) * Ebc;
            Ie = -(Ic + Ib);
            double sh = gadd * (Volt(v, d.A) - Volt(v, d.C));
            Ic += sh; Ie -= sh;
        }

        // чистая функция: невязка Кирхгофа в точке v, ничего не меняет
        double[] Residual(double[] v, double gadd)
        {
            var r = new double[M];
            for (int i = 0; i < nl.Devices.Count; i++)
            {
                var d = nl.Devices[i];
                int a = idx(d.A), b = idx(d.B);
                if (d.Kind == Kind.R)
                {
                    double cur = (Volt(v, d.A) - Volt(v, d.B)) / d.Value;
                    if (a >= 0) r[a] += cur;
                    if (b >= 0) r[b] -= cur;
                }
                else if (d.Kind == Kind.V)
                {
                    int m = branch[i];
                    if (a >= 0) r[a] += v[m];
                    if (b >= 0) r[b] -= v[m];
                    r[m] += Volt(v, d.A) - Volt(v, d.B) - d.Value;
                }
                else
                {
                    TrueBJT(d, v, gadd, out double Ic, out double Ib, out double Ie);
                    int c = idx(d.A), bb = idx(d.B), e = idx(d.C);
                    if (c >= 0) r[c] += Ic;
                    if (bb >= 0) r[bb] += Ib;
                    if (e >= 0) r[e] += Ie;
                }
            }
            return r;
        }

        // сборка линеаризованной системы в точке v (обновляет prev)
        void Stamp(double[] v, double gadd, double[,] J, double[] rhs)
        {
            double vcrit = Vt * Math.Log(Vt / (Is * Math.Sqrt(2)));
            for (int i = 0; i < nl.Devices.Count; i++)
            {
                var d = nl.Devices[i];
                int a = idx(d.A), b = idx(d.B);
                if (d.Kind == Kind.R) g2(J, a, b, 1 / d.Value);
                else if (d.Kind == Kind.V)
                {
                    int m = branch[i];
                    if (a >= 0) { J[a, m] += 1; J[m, a] += 1; }
                    if (b >= 0) { J[b, m] -= 1; J[m, b] -= 1; }
                    rhs[m] += d.Value;
                }
                else
                {
                    int c = idx(d.A), bb = idx(d.B), e = idx(d.C);
                    double aF = d.Value / (d.Value + 1);
                    double bR = 0.1 * d.Value;
                    double aR = Math.Max(0.5, bR / (bR + 1));
                    double vbe = Volt(v, d.B) - Volt(v, d.C), vbc = Volt(v, d.B) - Volt(v, d.A);
                    if (!prev.ContainsKey(i)) prev[i] = new double[] { vbe, vbc };
                    vbe = pnLim(vbe, prev[i][0], vcrit);
                    vbc = pnLim(vbc, prev[i][1], vcrit);
                    prev[i][0] = vbe; prev[i][1] = vbc;
                    double Ebe = Is * (Math.Exp(vbe / Vt) - 1);
                    double Ebc = Is * (Math.Exp(vbc / Vt) - 1);
                    double Ic = aF * Ebe - Ebc;
                    double Ib = (1 - aF) * Ebe + (1 - aR) * Ebc;
                    double Ie = -(Ic + Ib);
                    double Gee = Math.Max((Ebe + Is) / Vt, Gmin);
                    double Gcc = Math.Max((Ebc + Is) / Vt, Gmin);
                    double Gce = -Gee * aF, Gec = -Gcc * aR;
                    g2(J, bb, e, Gee);
                    g2(J, c, bb, Gcc);
                    cross(J, bb, c, e, Gce);
                    cross(J, bb, e, c, Gec);
                    g2(J, c, e, gadd);
                    if (bb >= 0) rhs[bb] -= Ib - (Gcc + Gec) * vbc - (Gee + Gce) * vbe;
                    if (c >= 0) rhs[c] -= Ic + Gcc * vbc + Gce * vbe;
                    if (e >= 0) rhs[e] -= Ie + Gee * vbe + Gec * vbc;
                }
            }
        }

        public double[] Solve()
        {
            var x = new double[M];
            for (int it = 0; it < 300; it++)
            {
                double gadd = it > 100 ? 1e-3 : 1e-6;
                var J = new double[M, M];
                var rhs = new double[M];
                Stamp(x, gadd, J, rhs);                    // линеаризуем в текущей x
                double norm = MaxAbs(Residual(x, gadd));
                if (norm < 1e-7) break;                    // сошлось
                double[] xn = Gauss(J, rhs);               // решили линейную систему
                double alpha = 1;
                while (alpha > 1e-3)                       // line search: пробуем шаг, не портя x
                {
                    var trial = new double[M];
                    for (int k = 0; k < M; k++) trial[k] = x[k] + alpha * (xn[k] - x[k]);
                    if (MaxAbs(Residual(trial, gadd)) < norm) { x = trial; break; }
                    alpha *= 0.5;
                }
                Iters = it;
                FinalNorm = norm;
            }
            return x;
        }
    }
}