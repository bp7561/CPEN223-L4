// Lab 4
// Student name: Brian Pham
// Student number: 24626509

using System;
using System.Collections.Generic;

Console.WriteLine("CPEN223 Lab 4");

//Testing: Write test cases that exercise all four methods you are to implement.
//TODO
Console.WriteLine("Default test cases:");
Dictionary<string, int> counts = GenomeAnalyzer.CountKMers("AAAA", 2);
Console.WriteLine($"Expected: AA -> 3, Actual: AA -> {counts["AA"]} ({counts.Count} entries)");

bool differ = GenomeAnalyzer.SamplesDiffer("AAAA", "TTTT", 2, 3);
Console.WriteLine($"Expected: True, Actual: {differ}");
Console.WriteLine("");

Console.WriteLine("k = 1 test case:");
Dictionary<string, int> countsk = GenomeAnalyzer.CountKMers("AACGT", 1); // A -> 2, C -> 1, G -> 1, T -> 1
Console.WriteLine($"Expected: A -> 2, C -> 1, G -> 1, T -> 1, Actual: A -> {countsk["A"]}, C -> {countsk["C"]}, G -> {countsk["G"]}, T -> {countsk["T"]} ({countsk.Count} entries)");
Console.WriteLine("");

Console.WriteLine("k equal to the full sequence length test case:");
Dictionary<string, int> countskf = GenomeAnalyzer.CountKMers("ACGT", 4); // ACGT -> 1
Console.WriteLine($"Expected: ACGT -> 1, Actual: ACGT -> {countskf["ACGT"]} ({countskf.Count} entries)");
Console.WriteLine("");

Console.WriteLine("A k-mer present in only one sample test case:");
Dictionary<string, int> changeone = GenomeAnalyzer.CompareProfiles("AAAA", "TTTT", 2); // AA -> -3, TT -> +3
Console.WriteLine($"Expected: AA -> -3, TT -> +3, Actual: AA -> {changeone["AA"]}, TT -> {changeone["TT"]} ({changeone.Count} entries)");
Console.WriteLine("");

Console.WriteLine("A tie for the largest change test case:");
List<string> resulttie = GenomeAnalyzer.MostChangedKMers("AAAA", "TTTT", 2); // returns AA and TT (order does not matter)
Console.WriteLine($"Expected: 2, Actual: ({resulttie.Count} entries)");
Console.WriteLine("");

Console.WriteLine("Identical samples test cases:");
Dictionary<string, int> changeidentical = GenomeAnalyzer.CompareProfiles("ACGTAC", "ACGTAC", 3); // empty dictionary
Console.WriteLine($"Expected: 0, Actual: ({changeidentical.Count} entries)");	
List<string> resultidentical = GenomeAnalyzer.MostChangedKMers("ACGTAC", "ACGTAC", 3); // empty list
Console.WriteLine($"Expected: 0, Actual: ({resultidentical.Count} entries)");
Console.WriteLine("");
//end Testing code

//Do not change the program skeleton: keep the class name, method names,
//parameters, and return types exactly as given.
//Do not use LINQ, and do not use Console inside the GenomeAnalyzer methods.

public static class GenomeAnalyzer
{
	// helper function to validate sequences
	private static void CheckSequence(string sequence)
    {
		// sequence is empty
        if (sequence == null)
        {
            throw new ArgumentException("ERROR: sequence is null");
        }
		
		// iterate through each character in sequence 
        foreach (char c in sequence)
        {
			// ACGT only accepted
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            {
                throw new ArgumentException("ERROR: sequence is invalid");
            }
        }
    }
	
    public static Dictionary<string, int> CountKMers(string sequence, int k)
    {
		// validation
		CheckSequence(sequence);
		
        if (k <= 0)
		{
			throw new ArgumentException("ERROR: k is lesser than or equal to 0");
		}
		if (k > sequence.Length)
		{
			throw new ArgumentException("ERROR: k is greater than sequence length");
		}
		
		// new empty dictionary to hold kmer count
		Dictionary<string, int> count = new Dictionary<string, int>();
		
		// sliding window to iterate through sequence to check for kmers
		for (int i = 0; i <= sequence.Length - k; i++) 
		{
			// extract string of length k from sequence
			string kmer = sequence.Substring(i, k);
			
			// checks if dictionary has the kmer
			if (count.TryGetValue(kmer, out int currCount))
			{
				// update kmer count
				count[kmer] = currCount + 1;	
			}
			// else note new kmer
			else
			{
				count[kmer] = 1;
			}
		}
		// return dictionary
		return count;
		
    }

    public static Dictionary<string, int> CompareProfiles(string reference, string sample, int k)
    {
        // validation
		CheckSequence(reference);
		CheckSequence(sample);
		
        if (k <= 0)
		{
			throw new ArgumentException("ERROR: k is lesser than or equal to 0");
		}
		if (k > reference.Length || k > sample.Length)
		{
			throw new ArgumentException("ERROR: k is greater than reference or sample length");
		}
		
		// call CountKmers for kmer count in reference and sample sequence
		Dictionary<string, int> refCount = CountKMers(reference, k);
		Dictionary<string, int> sampleCount = CountKMers(sample, k);
		
		// new empty dictionary to hold difference between profiles count
		Dictionary<string, int> diff = new Dictionary<string, int>();
		
		// iterate through reference sequence
		foreach (KeyValuePair<string, int> pair in refCount)
		{
			// seperate kmer and its count
			string kmer = pair.Key;
			int refCounter = pair.Value;
			
			// checks if kmer exists in sample sequence, if none, then 0
			sampleCount.TryGetValue(kmer, out int sampleCounter);
			
			// calculate difference
			int change = sampleCounter - refCounter;
			
			// if there is a difference, note difference
			if (change != 0)
			{
				diff[kmer] = change;
			}
		}
		
		// iterate through sample sequence
		foreach (KeyValuePair<string, int> pair in sampleCount)
		{
			// seperate kmer and its count
			string kmer = pair.Key;
			int sampleCounter = pair.Value;
			
			// checks if kmer is not noted in reference
			if (!refCount.ContainsKey(kmer))
			{
				// calculate difference without reference counter
				int change = sampleCounter;
				
				// note difference
				if (change != 0)
				{
					diff[kmer] = change;
				}
			}
		}
		
		// return dictionary
		return diff;
    }

    public static List<string> MostChangedKMers(string reference, string sample, int k)
    {
		// call CompareProfiles for differences in sequences
        Dictionary<string, int> change = CompareProfiles(reference, sample, k);
		// new empty list to hold most changed kmers
		List<string> result = new List<string>(); 
		
		// initialize max changed noted to 0
		int maxChanged = 0;
		
		// iterate through dictionary with kmer differences
		foreach (KeyValuePair<string, int> pair in change)
		{
			// deal with positive differences only
			int absChanged = Math.Abs(pair.Value);
			
			// if kmer difference is greater than current change, replace
			if (absChanged > maxChanged)
			{
				maxChanged = absChanged;
				// clear current list and add new kmer
				result.Clear();
				result.Add(pair.Key);
			}
			// else if same, add to list
			else if (absChanged == maxChanged && maxChanged > 0)
			{
				result.Add(pair.Key);
			}
		}
		// return list
		return result;
    }

    public static bool SamplesDiffer(string reference, string sample, int k, int threshold)
    {
        // validation
        if (threshold <= 0)
        {
            throw new ArgumentException("ERROR: threshold is lesser than or equal to 0.");
        }

        // call CompareProfiles for differences in sequences
        Dictionary<string, int> change = CompareProfiles(reference, sample, k);

        // iterate through dictionary with kmer differences
        foreach (KeyValuePair<string, int> pair in change)
        {
			// if differences exceed threshold
            if (Math.Abs(pair.Value) >= threshold)
            {
                return true;
            }
        }
		// else 
        return false;
    }
}
