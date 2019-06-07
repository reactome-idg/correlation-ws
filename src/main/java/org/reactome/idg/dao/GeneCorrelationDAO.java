package org.reactome.idg.dao;

import java.util.Map;

import org.reactome.idg.model.Provenance;

public interface GeneCorrelationDAO
{
	public void setCurrentProvenance(Provenance p);
	
	public Provenance getCurrentProvenance();
	
	public void addGenePair(String gene1, String gene2, double correlationValue);
	
	public Provenance addProvenance(Provenance p);
	
	public Map<Provenance, Double> getCorrelation(String gene1, String gene2);
	
	public Double getCorrelation(String gene1, String gene2, Provenance prov);
	
	public void setBatchSize(int batchSize);
	
	public int getBatchSize();
	
	public int getNumTxOps();
	
	public void loadGenePairsFromDataFile(String pathToFile);
}
