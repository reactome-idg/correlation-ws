package org.reactome.idg;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.reactome.idg.loader.HarmonizomeLoader;

/**
 * Contains data loaded from various Harmonizome data sources.
 * @author sshorser
 *
 */
public class HarmonizomeDataRepository implements DataRepository
{

	// A mapping of a Provenance to a loader.
	private Map<Provenance, HarmonizomeLoader> dataLoaders;
	// A mapping of Provenance to actual data.
	private Map<Provenance, TreeMap<String,List<Double>>> dataByGeneSymbol = new TreeMap<>();
	private Map<Provenance, TreeMap<String,List<Double>>> dataByUniprotAccession = new TreeMap<>();;

	
	/**
	 * Private constructor. Use the factory method <code>createDataRepository</code>.
	 * @param loaders
	 */
	private HarmonizomeDataRepository(Map<Provenance, HarmonizomeLoader> loaders)
	{
		this.dataLoaders = loaders;
	}
	
	/**
	 * Creates a new DataLoader for Harmonizome.
	 * @param loaders - a map of loaders, keyed by their provenance.
	 * @return A new HarmonizomeDataRepository.
	 */
	public static HarmonizomeDataRepository createDataRepository(Map<Provenance, HarmonizomeLoader> loaders)
	{
		HarmonizomeDataRepository repository = new HarmonizomeDataRepository(loaders);
		
		return repository;
	}
	
	public void executeDataLoaders()
	{
		for (Provenance p : this.dataLoaders.keySet())
		{
			HarmonizomeLoader loader = this.dataLoaders.get(p);
			try
			{
				long numRecords = loader.loadGeneAssociationData();
				this.dataByGeneSymbol.put(p, loader.getGeneSymbolDataset());
				this.dataByUniprotAccession.put(p, loader.getUniprotAccessionDataset());
			}
			catch (FileNotFoundException e)
			{
				e.printStackTrace();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Returns the correlations between the two genes, mapped by their provenance.
	 * @param gene1
	 * @param gene2
	 * @return A map, keyed by Provenance.
	 */
	@Override
	public Map<Provenance, List<Double>> getGeneCorrelation(String gene1, String gene2, LookupKeyType keyType)
	{
		String key = DataRepository.generateKey(gene1, gene2);
		Map<Provenance, List<Double>> correlations = new TreeMap<>();
		if (keyType == LookupKeyType.GENE_SYMBOL)
		{
			for (Provenance p : dataByGeneSymbol.keySet())
			{
				if (this.dataByGeneSymbol.get(p).containsKey(key))
				{
					correlations.put(p, this.dataByGeneSymbol.get(p).get(key));
				}
			}
		}
		else
		{
			for (Provenance p : dataByUniprotAccession.keySet())
			{
				if (this.dataByUniprotAccession.get(p).containsKey(key))
				{
					correlations.put(p, this.dataByUniprotAccession.get(p).get(key));
				}
			}
		}
		return correlations;
	}
	
	/**
	 * Returns the correlation between two genes for a specific provenance. A list is returned since it's possible that GeneSymbol:UniProt mappings
	 * are not 1:1, which means there could be gene-pairs with more than 1 correlation.
	 * @param prov
	 * @param gene1
	 * @param gene2
	 * @return The correlation value.
	 */
	@Override
	public List<Double> getGeneCorrelation(Provenance prov, String gene1, String gene2, LookupKeyType keyType)
	{
		List<Double> correlation;
		String key = DataRepository.generateKey(gene1, gene2);
		if (keyType == LookupKeyType.GENE_SYMBOL)
		{
			correlation = this.dataByGeneSymbol.get(prov).get(key);
		}
		else
		{
			correlation = this.dataByUniprotAccession.get(prov).get(key);
		}
		return correlation;
	}

	@Override
	public Map<Provenance, List<Double>> getGeneCorrelation(String gene1, String gene2)
	{
		// If lookup-key type is not known, we'll just have to try both Gene Symbol and if that fails then try UniProt...
		Map<Provenance, List<Double>> correlations = this.getGeneCorrelation(gene1, gene2, LookupKeyType.GENE_SYMBOL);
		if (correlations == null || correlations.isEmpty())
		{
			correlations = this.getGeneCorrelation(gene1, gene2, LookupKeyType.UNIPROT_ACCESSION);
		}
		return correlations;
	}

	@Override
	public List<Double> getGeneCorrelation(Provenance prov, String gene1, String gene2)
	{
		// If lookup-key type is not known, we'll just have to try both Gene Symbol and if that fails then try UniProt...
		List<Double> correlation = this.getGeneCorrelation(prov, gene1, gene2, LookupKeyType.GENE_SYMBOL);
		if (correlation == null)
		{
			correlation = this.getGeneCorrelation(prov, gene1, gene2, LookupKeyType.UNIPROT_ACCESSION);
		}
		return correlation;
	}
	
}
