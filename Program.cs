using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace BJTgateLab
{
    enum Kind { V, R, Q }                 // тип элемента

    class Device                          // один элемент схемы
    {
        public string Name = "";
        public Kind Kind;
        public int A, B, C;               // узлы: R/V -> A,B; Q -> A=C(колл), B=база, C=эмит
        public double Value;              // R: ом, V: вольт, Q: beta
    }

    class Netlist                         // вся схема
    {
        public List<Device> Devices = new List<Device>();
        public List<string> NodeNames = new List<string>();
        public Dictionary<string, int> NodeIndex = new Dictionary<string, int>();

        public Netlist() { NodeNames.Add("0"); NodeIndex["0"] = 0; }  // земля = узел 0

        public int Node(string name)      // имя узла -> номер (новые добавляются)
        {
            if (!NodeIndex.ContainsKey(name))
            {
                NodeIndex[name] = NodeNames.Count;
                NodeNames.Add(name);
            }
            return NodeIndex[name];
        }
    }

    class Program
    {
        static Netlist Parse(string[] lines)
        {
            var nl = new Netlist();
            foreach (var raw in lines)
            {
                string line = raw.Trim();
                if (line.Length == 0 || line[0] == '#' || line[0] == '*') continue;
                // твой же приём из Timus:
                string[] t = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                var d = new Device { Name = t[1] };
                switch (t[0].ToUpper())
                {
                    case "V":
                        d.Kind = Kind.V;
                        d.A = nl.Node(t[2]); d.B = nl.Node(t[3]);
                        d.Value = double.Parse(t[4], CultureInfo.InvariantCulture);
                        break;
                    case "R":
                        d.Kind = Kind.R;
                        d.A = nl.Node(t[2]); d.B = nl.Node(t[3]);
                        d.Value = double.Parse(t[4], CultureInfo.InvariantCulture);
                        break;
                    case "Q":
                        d.Kind = Kind.Q;
                        d.A = nl.Node(t[2]); d.B = nl.Node(t[3]); d.C = nl.Node(t[4]);
                        d.Value = double.Parse(t[5], CultureInfo.InvariantCulture);
                        break;
                    default:
                        throw new Exception("Неизвестный элемент: " + t[0]);
                }
                nl.Devices.Add(d);
            }
            return nl;
        }

        static void Main(string[] args)
        {
            string path = args.Length > 0 ? args[0] : "nand.txt";
            Netlist nl = Parse(File.ReadAllLines(path));

            Console.WriteLine("Узлы:");
            for (int i = 0; i < nl.NodeNames.Count; i++)
                Console.WriteLine($"  {i} = {nl.NodeNames[i]}");

            Console.WriteLine("Элементы:");
            foreach (var d in nl.Devices)
            {
                string nodes = d.Kind == Kind.Q ? $"{d.A},{d.B},{d.C}" : $"{d.A},{d.B}";
                Console.WriteLine($"  {d.Name} [{d.Kind}] узлы=({nodes}) значение={d.Value}");
            }
        }
    }
}