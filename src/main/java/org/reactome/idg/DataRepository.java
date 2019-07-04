package org.reactome.idg;

import java.util.List;
import java.util.Map;
import org.reactome.idg.model.Provenance;

/**
 * Interface for a Reactome-IDG data repository. A repository can contain a collection of datasets, identified by their provenance.
 * The idea is to have one DataRepository for each major data source/distinct Provenance, i.e. one DataRepository for Harmonizome, one for ARCHS4, etc...
 * @author sshorser
 *
 */
// I now realize that the amount of data is probably WAY too much to store in memory, so this interface and its implementations may change significantly or
// go away entirely.
public interface DataRepository
{
	// We may not be storing data keyed by UniProt accession so this enum may eventually be removed.
	public enum LookupKeyType
	{
		GENE_SYMBOL, UNIPROT_ACCESSION;
	}
	
	public Map<Provenance, List<Double>> getGeneCorrelation(String gene1, String gene2, LookupKeyType keyType);
	
	public List<Double> getGeneCorrelation(Provenance prov, String gene1, String gene2, LookupKeyType keyType);

	public Map<Provenance, List<Double>> getGeneCorrelation(String gene1, String gene2);
	
	public List<Double> getGeneCorrelation(Provenance prov, String gene1, String gene2);
	
	static String generateKey(String symbol1, String symbol2)
	{
		return symbol1.compareTo(symbol2) <= 0
				? symbol1 + "|" + symbol2
				: symbol2 + "|" + symbol1;
	}
}
