package org.reactome.idg.loader;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.reactome.idg.DataRepository;

// TODO: Harmonizome data should be loaded directly into MySQL.

/**
 * Loads data from a Harmonizome Gene Association file. This loads data into in-memory
 * data structures.
 * 
 * This class may become obsolete if we don't need the Uniprot Accessions that are
 * in the Harmonizome file.
 * @author sshorser
 *
 */
public class HarmonizomeLoader
{
//	private Provenance dataProvenance;
	private String filePath;
	private TreeMap<String, List<Double>> correlationValuesByGeneSymbol;
	private TreeMap<String, List<Double>> correlationValuesByUniprotAcc;
	
	public HarmonizomeLoader(String filePath/* , Provenance provenance */)
	{
//		this.dataProvenance = provenance;
		this.filePath = filePath;
	}
	
	public TreeMap<String, List<Double>> getGeneSymbolDataset()
	{
		return this.correlationValuesByGeneSymbol;
	}
	
	public TreeMap<String, List<Double>> getUniprotAccessionDataset()
	{
		return this.correlationValuesByUniprotAcc;
	}
	
	/**
	 * Loads the Gene Association data, returns it in a structure.
	 * @return the number of distinct gene-pair correlations loaded.
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public long loadGeneAssociationData() throws FileNotFoundException, IOException
	{
		// Harmonizome files contain 3 lines of headers: the first has gene symbols, the second has
		// corresponding Uniprot accessions, the third is the NCBI gene ID (I assume an internal Harminozome identifier).
		//
		// I think FileInputStream with Scanner will allow us to read the file while only keeping one line
		// in memory at a time, which is important if the files are large (ARCHS4 files are large: ~13GB as CSV!) 
		try(FileInputStream inStream = new FileInputStream(this.filePath);
			Scanner scanner = new Scanner(inStream))
		{
			// We need to store the Gene Symbol and UniProt Accession.
			/*
			 * The file looks a bit like this:

#         #             GeneSym       PRKAR1A     PRKAR1B     AMICA1    ADAM28    PTPRA...
#         #             UniProtACC    P10644      P31321      Q86YT9    Q9UKQ2    P18433...
GeneSym   UniProtACC    GeneID/GeneID 5573        5575        120425    10863     5786
PRKAR1A   P10644        5573          1.000000    1.000000    1.000...
PRKAR1B   P31321        5575          1.000000    1.000000    1.000...
AMICA1    Q86YT9        120425        1.000000    1.000000    1.000...
ADAM28    Q9UKQ2        10863         1.000000    1.000000    1.000...
PTPRA     P18433        5786          0.707107    0.707107    0.707...
...

			 */
			if (scanner.hasNext())
			{
				// Header 1 is Gene Symbols. 
				String geneSymbolHeader = scanner.nextLine();
				String[] parts = geneSymbolHeader.split("\t");
				String[] allGeneSymbols = new String[parts.length - 3];
				for (int i = 3; i < parts.length; i++)
				{
					allGeneSymbols[i - 3] = parts[i];
				}
				System.out.println(allGeneSymbols.length + " gene symbols in the header.");
				System.out.println( (allGeneSymbols.length * (allGeneSymbols.length + 1)/2) + " possible gene-pair correlations.");
				// Line 2 is UniProt accessions
				String uniprotHeader = scanner.nextLine();
				parts = uniprotHeader.split("\t");
				String[] allUniprotAccessions = new String[parts.length - 3];
				for (int i = 3; i < parts.length; i++)
				{
					allUniprotAccessions[i - 3] = parts[i];
				}
				
				// Line 3 is NCBI Gene IDs - not used in this program.
				scanner.nextLine();
				
				LocalDateTime startTime = LocalDateTime.now();
				
				// Now that we know all the gene symbols we could get, we can start building the mapping.
				// The internal structure will be a map of gene-pair to correlation value. 
				// Keys will be of the form <GeneSymbol1>-<GeneSymbol2> - with GeneSymbol1 and GeneSymbol always in alphabetical sort-order.
				//
				// we don't need to start reading each line at the leftmost column, since this should be a triangular matrix.
				// Each line down, we can start reading one position further to the right. The point of this is to *avoid* processing every gene-pair twice.
				int lineStartOffset = 3;
				correlationValuesByGeneSymbol = new TreeMap<>();
				correlationValuesByUniprotAcc = new TreeMap<>();
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					parts = line.split("\t");
					String currentGeneSymbol = parts[0];
					String currentUniprotAcc = parts[1];
					
					//TODO: process the line in two threads: one for gene symbol, one for uniprot acc.
					// Eventually, some of this will probably need to be rewritten for better performance. When I added
					// code to process UniProt data, it almost doubled the running time (for larger - 55MB - files) to 10 seconds.
					// 10 seconds isn't TOO bad, but there are much larger files out there...
					for (int i = lineStartOffset; i < parts.length; i++)
					{
						int lookupIndex = i - 3;
						String correlationValue = parts[i];
						// Get the other gene symbol for this position in the line.
						String otherGeneSymbol = allGeneSymbols[lookupIndex];
						String key = DataRepository.generateKey(currentGeneSymbol, otherGeneSymbol);
						final Double correlationAsDouble = Double.valueOf(correlationValue);
						if (!correlationValuesByGeneSymbol.containsKey(key))
						{
							correlationValuesByGeneSymbol.put(key, new ArrayList<>( Arrays.asList(correlationAsDouble)) );
						}
						// If by some chance this key has already been seen (possibly because of a duplicated identifier), append the correlation value.
						else 
						{
							List<Double> values = correlationValuesByGeneSymbol.get(key);
							values.add(correlationAsDouble);
							System.out.println("Whoa! The GeneSymbol key " + key + " has multiple values: " + String.join(",", values.stream().map( v -> v.toString()).collect(Collectors.toList()) ) );
							correlationValuesByGeneSymbol.put(key, values);
						}
						
						// Get the other gene symbol for this position in the line.
						String otherUniprotAcc = allUniprotAccessions[lookupIndex];
						String uniprotKey = DataRepository.generateKey(currentUniprotAcc, otherUniprotAcc);
						if (!correlationValuesByUniprotAcc.containsKey(uniprotKey))
						{
							correlationValuesByUniprotAcc.put(uniprotKey, new ArrayList<>( Arrays.asList(correlationAsDouble)) );
						}
						// If by some chance this key has already been seen (possibly because of a duplicated identifier), append the correlation value.
						else 
						{
							List<Double> values = correlationValuesByUniprotAcc.get(uniprotKey);
							values.add(correlationAsDouble);
							System.out.println("Whoa! The UniProt key " + uniprotKey + " has multiple values: " + String.join(",", values.stream().map( v -> v.toString()).collect(Collectors.toList()) ) );
							correlationValuesByUniprotAcc.put(uniprotKey, values);
						}
					}
					lineStartOffset++;
				}
				LocalDateTime endTime = LocalDateTime.now();
				// Now that we're done, print out the elapsed time and number of keys loaded.
				System.out.println(correlationValuesByGeneSymbol.keySet().size() + " gene-pair keys (Gene Symbol) have been loaded.");
//				System.out.println(String.join(",", correlationValuesByGeneSymbol.keySet()));
				System.out.println(correlationValuesByUniprotAcc.keySet().size() + " gene-pair keys (Uniprot Accession) have been loaded.");
//				System.out.println(String.join(",", correlationValuesByUniprotAcc.keySet()));
				System.out.println(Duration.between(startTime, endTime).toString() + " time spent loading the data.");
			}
		}
		return correlationValuesByGeneSymbol.keySet().size();
	}



//	public Provenance getDataProvenance()
//	{
//		return dataProvenance;
//	}
}
