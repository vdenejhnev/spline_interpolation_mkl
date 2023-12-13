using System;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;

public struct DataItem
{
    public double X;
    public double[] Y;

    public DataItem(double x, double y1, double y2)
    {
        X = x;
        Y = new double[2] { y1, y2 };
    }

    public string ToLongString(string format)
    {
        return $"X: {X.ToString(format)}, Y: [{Y[0].ToString(format)} , {Y[1].ToString(format)}]";
    }

    public override string ToString()
    {
        return $"X: {X}, Y: [{Y[0]} , {Y[1]}]";
    }
}

public abstract class V2Data : IEnumerable<DataItem>
{
    public string Key { get; set; }
    public DateTime Time { get; set; }

    public V2Data(string key, DateTime time)
    {
        Key = key;
        Time = time;
    }

    public abstract IEnumerator<DataItem> GetEnumerator();

    System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public abstract double MinField { get; }

    public abstract bool IsNull { get; }

    public abstract string ToLongString(string format);

    public override string ToString()
    {
        return $"Key: {Key}, Time: {Time.ToString("yyyy-MM-dd HH:mm:ss")}";
    }
}

public class V2DataList : V2Data
{
    public delegate double[] FDI(double x);

    public List<DataItem> DataItems { get; set; }

    public V2DataList(string key = "default", DateTime date = new System.DateTime()) : base(key, date)
    {
        DataItems = new List<DataItem>();
    }

    public V2DataList(string key, DateTime date, double[] x, FDI F) : base(key, date)
    {
        DataItems = new List<DataItem>();
        for (int i = 0; i < x.Length; i++)
        {
            double[] y = F(x[i]);
            DataItem dataItem = new DataItem(x[i], y[0], y[1]);
            if (!DataItems.Any(item => item.X == x[i]))
            {
                DataItems.Add(dataItem);
            }
        }
    }

    public override double MinField
    {
        get
        {
            double min_field = DataItems[0].Y[0];
            foreach (DataItem data_item in DataItems)
            {
                if (data_item.Y[0] < min_field)
                {
                    min_field = data_item.Y[0];
                }

                if (data_item.Y[1] < min_field)
                {
                    min_field = data_item.Y[0];
                }
            }
            return min_field;
        }
    }

    public override bool IsNull
    {
        get
        {
            foreach (var item in DataItems)
            {
                if (item.Y[0] == 0 && item.Y[1] == 0)
                {
                    return true;
                }
            }
            return false;
        }
    }

    public V2DataArray DataArray
    {
        get
        {
            V2DataArray dataArray = new V2DataArray(this.Key, this.Time);
            return dataArray;
        }
    }

    public override string ToString()
    {
        return $"V2DataList, {base.ToString()}, Number of DataItems: {DataItems.Count}";
    }

    public override string ToLongString(string format)
    {
        string result = "";
        result += $"{ToString()}\n";
        result += "DataItems:\n";
        foreach (var item in DataItems)
        {
            result += $"{item.ToLongString(format)}\n";
        }
        return result;
    }

    public void Add(DataItem dataItem)
    {
        if (!DataItems.Any(item => item.X == dataItem.X))
        {
            DataItems.Add(dataItem);
        }
    }

    public override IEnumerator<DataItem> GetEnumerator()
    {
        return DataItems.GetEnumerator();
    }
}

public class V2DataArray : V2Data
{
    public delegate double FValues(double x, int index);

    public double[] X { get; set; }

    public double[,] Y { get; set; }

    public V2DataArray(string key = "default", DateTime date = new System.DateTime()) : base(key, date)
    {
        X = new double[0];
        Y = new double[0, 0];
    }

    public V2DataArray(string key, DateTime date, double[] x, FValues F) : base(key, date)
    {
        X = x;
        Y = new double[x.Length, 2];

        for (int i = 0; i < x.Length; i++)
        {
            Y[i, 0] = F(x[i], 0);
            Y[i, 1] = F(x[i], 1);
        }
    }

    public V2DataArray(string key, DateTime date, int nX, double xL, double xR, FValues F) : base(key, date)
    {
        X = new double[nX];
        Y = new double[nX, 2];

        double dx = (xR - xL) / (nX - 1);

        for (int i = 0; i < nX; i++)
        {
            X[i] = xL + i * dx;
            Y[i, 0] = F(X[i], 0);
            Y[i, 1] = F(X[i], 1);
        }
    }

    public double[] this[int index]
    {
        get
        {
            if (index == 0 && index == 1)
            {
                return Enumerable.Range(0, X.Length).Select(i => Y[i, index]).ToArray();
            }
            else
            {
                throw new ArgumentOutOfRangeException(nameof(index));
            }
        }
    }

    public override double MinField => Y.Cast<double>().Min(Math.Abs);

    public override bool IsNull
    {
        get
        {
            for (int i = 0; i < X.Length; i++)
            {
                if (Y[i, 0] == 0 && Y[i, 1] == 0)
                {
                    return true;
                }
            }

            return false;
        }
    }

    public static explicit operator V2DataList(V2DataArray source)
    {
        V2DataList result = new V2DataList(source.Key, source.Time);

        for (int i = 0; i < source.X.Length; i++)
        {
            var dataItem = new DataItem(source.X[i], source.Y[i, 0], source.Y[i, 1]);
            result.Add(dataItem);
        }

        return result;
    }

    public override string ToString()
    {
        return $"V2DataArray, {base.ToString()}";
    }

    public override string ToLongString(string format)
    {
        string result = "";
        result += $"{ToString()}\n";
        result += "DataItems:\n";
        for (int i = 0; i < X.Length; i++)
        {
            result += $"X = {X[i].ToString(format)}, Y = ({Y[i, 0].ToString(format)}, {Y[i, 1].ToString(format)})\n";
        }

        return result;
    }

    public override IEnumerator<DataItem> GetEnumerator()
    {
        for (int i = 0; i < X.Length; i++)
        {
            var dataItem = new DataItem(X[i], Y[i, 0], Y[i, 1]);
            yield return dataItem;
        }
    }

    public static bool Save(string filename, V2DataArray dataArray)
    {
        StreamWriter writer = new StreamWriter(File.Open(filename, FileMode.Create));

        try
        {
            writer.WriteLine(dataArray.Key);
            writer.WriteLine(dataArray.Time);
            writer.WriteLine(dataArray.X.Length);
            writer.WriteLine();

            for (int i = 0; i < dataArray.X.Length; i++)
            {
                writer.WriteLine(dataArray.X[i].ToString());
                writer.WriteLine(dataArray.Y[i, 0].ToString());
                writer.WriteLine(dataArray.Y[i, 1].ToString());
                writer.WriteLine();
            }

            Console.WriteLine($"Success saving data to file '{filename}'");
            return true;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Error saving data to file '{filename}': {ex.Message}");
            return false;
        }
        finally
        {
            writer.Close();
        }
    }

    public static bool Load(string filename, ref V2DataArray dataArray)
    {
        StreamReader reader = new StreamReader(filename);
        string state = "S";
        int item = 0;

        try
        {
            while (!reader.EndOfStream)
            {
                switch (state)
                {
                    case "S":
                        dataArray.Key = reader.ReadLine();
                        dataArray.Time = Convert.ToDateTime(reader.ReadLine());
                        int length = Convert.ToInt32(reader.ReadLine());
                        dataArray.X = new double[length];
                        dataArray.Y = new double[length, 2];
                        reader.ReadLine();
                        state = "X";
                        break;
                    case "X":
                        dataArray.X[item] = Convert.ToDouble(reader.ReadLine());
                        state = "Y0";
                        break;
                    case "Y0":
                        dataArray.Y[item, 0] = Convert.ToDouble(reader.ReadLine());
                        state = "Y1";
                        break;
                    case "Y1":
                        dataArray.Y[item, 1] = Convert.ToDouble(reader.ReadLine());
                        state = "N";
                        break;
                    case "N":
                        reader.ReadLine();
                        item++;
                        state = "X";
                        break;
                }
            }

            Console.WriteLine($"Success loading data from file '{filename}'");

            return true;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Error loading data from file '{filename}': {ex.Message}");
            return false;
        }
        finally
        {
            reader.Close();
        }
    }
}

public class V2MainCollection : ObservableCollection<V2Data>
{

    public V2MainCollection()
    {
        List<V2Data> collectionList = new List<V2Data>();
    }

    public V2MainCollection(int nV2DataArray, int nV2DataList) : base()
    {
        for (int i = 0; i < nV2DataArray; i++)
        {
            Add(new V2DataArray());
        }

        for (int i = 0; i < nV2DataList; i++)
        {
            Add(new V2DataList());
        }
    }

    public new IEnumerator<V2Data> GetEnumerator() => base.Items.GetEnumerator();

    public List<bool> IsNullList
    {
        get
        {
            return base.Items.Select(i => i.IsNull).ToList();
        }
    }

    public bool Contains(string key)
    {
        foreach (var item in base.Items)
        {
            if (item.Key == key)
            {
                return true;
            }
        }

        return false;
    }

    public new bool Add(V2Data v2Data)
    {
        if (!Contains(v2Data.Key))
        {
            base.Add(v2Data);
            return true;
        };

        return false;
    }

    public string ToLongString(string format)
    {
        string output = "V2MainCollection:\n";
        int i = 1;
        foreach (var elem in base.Items)
        {
            output += $"Item_{i}: " + elem.ToLongString(format) + "\n";
            i += 1;
        }
        return output;
    }

    public override string ToString()
    {
        string output = "V2MainCollection:\n";
        int i = 1;
        foreach (var elem in base.Items)
        {
            output += $"Item_{i}: " + elem.ToString() + "\n";
            i += 1;
        }
        return output;
    }

    public int MaxResultsCount
    {
        get
        {
            if (base.Items.Count == 0) return -1;

            var resultsCounts = base.Items.Select(data =>
            {
                return data.Select(item => Math.Abs(Math.Sqrt(item.Y[0] * item.Y[0] + item.Y[1] * item.Y[1]))).Where(item => item == 0).Count();
            });

            return resultsCounts.Max();
        }
    }

    public DataItem? DataItemMax
    {
        get
        {
            if (base.Items.Count == 0) return null;

            return base.Items.Select(data =>
            {
                return data.Select(item => item).OrderByDescending(item => Math.Abs(item.Y[0]) + Math.Abs(item.Y[1])).FirstOrDefault();
            }).OrderByDescending(item => Math.Abs(item.Y[0]) + Math.Abs(item.Y[1])).FirstOrDefault();
        }
    }

    public IEnumerable<double> UniqueCoordinates
    {
        get
        {
            if (base.Items.Count == 0) return null;

            return (from item in base.Items from item1 in item select item1.X).Where(item2 => (from item in base.Items from item1 in item select item1.X).Count(x => x == item2) == 1);
        }
    }
}


public struct SplineDataItem
{
    public double X { get; }
    public double Y { get; }
    public double ApproximatedValue { get; set; }

    public SplineDataItem(double x, double y, double approximatedValue)
    {
        X = x;
        Y = y;
        ApproximatedValue = approximatedValue;
    }

    public override string ToString()
    {
        return ToString("F2");
    }

    public string ToString(string format)
    {
        return $"({X.ToString(format)}, {Y.ToString(format)}) -> {ApproximatedValue.ToString(format)}";
    }
}


class SplineData
{
    public V2DataArray DataArray { get; set; }
    public int NumberOfMeshNodes { get; set; }
    public double[] CubicSplineValues { get; set; }
    public int MaximumNumberOfIterations { get; set; }
    public int StoppingReason { get; set; }
    public double MinimumNonZeros { get; set; }
    public List<SplineDataItem> Results { get; set; } 

    public SplineData(V2DataArray dataArray, int numberOfMeshNodes, int maximumNumberOfIterations)
    {
        DataArray = dataArray;
        NumberOfMeshNodes = numberOfMeshNodes;
        MaximumNumberOfIterations = maximumNumberOfIterations;
        Results = new List<SplineDataItem>();
    }

    public void ApproximateSpline()
    {
        double[] X = {DataArray.X[0], DataArray.X[DataArray.X.Length - 1]};
        double[] Y = new double[DataArray.X.Length];
        double[] sites = { DataArray.X[0], DataArray.X[DataArray.X.Length - 1] };
        int[] dorder = { 1 };
        double[] result = new double[NumberOfMeshNodes];
        int ret = 0;
        double[] scoeff = new double[4 * (DataArray.X.Length - 1)];

        for (int i = 0; i < Y.Length; i++)
        {
            Y[i] = DataArray.Y[i, 0];
        }

        SplineInterpolation(DataArray.X.Length, DataArray.X.Length, X, Y, scoeff, NumberOfMeshNodes, sites, 1, dorder, result, ref ret);

        if (ret == -1)
        {
            return;
        }

        Console.WriteLine("Results:");

        double rsn = 0;

        for (int i = 0; i < NumberOfMeshNodes; i++)
        {
            double site = X[0] + (i) * ((X[1] - X[0]) / (NumberOfMeshNodes - 1));
            double y = site * site;
            Results.Add(new SplineDataItem(site, y, result[i]));
            Console.WriteLine($"{site}:\t{y} | {result[i]}");
        }

        Console.WriteLine($"rsn: {rsn}");
    }

    public string ToLongString(string format)
    {
        StringBuilder result = new StringBuilder();
        result.AppendLine(DataArray.ToLongString(format));
        foreach (var resultItem in Results)
        {
            result.AppendLine(resultItem.ToString("f5"));
        }
        /*
        result.AppendLine($"Number of mesh nodes: {NumberOfMeshNodes}");
        result.AppendLine($"Minimum non-zero value: {MinimumNonZeros}");
        result.AppendLine($"Stopping reason: {StoppingReason}");
        result.AppendLine($"Number of iterations: {MaximumNumberOfIterations}");
        result.AppendLine("Spline data:");
        */
        return result.ToString();
    }

    public void Save(string filename, string format)
    {
        using (StreamWriter writer = new StreamWriter(filename))
        {
            writer.Write(ToLongString(format));
        }
    }

    [DllImport("..\\..\\..\\..\\x64\\Debug\\dll3.dll", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SplineInterpolation(int nx, int ny, double[] x, double[] y, double[] scoeff, int nsite, double[] site, int ndorder, int[] dorder, double[] result, ref int ret);
}


class Program
{
    static void Main(string[] args)
    {
        V2DataArray dataArray = new V2DataArray("DataArray", DateTime.Now, 10, 0, 100, (x, index) => x * x);
        Console.WriteLine(dataArray.ToLongString("f2"));
        SplineData splineData = new SplineData(dataArray, 6, 100);
        splineData.ApproximateSpline();
        splineData.Save("spline_data.txt", "f5");
    }
}