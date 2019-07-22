package org.reactome.idg.loader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;

public class Archs4ExpressionDataLoader
{
	private static final Logger logger = LogManager.getLogger();
	Archs4ExpressionMetadata metadata = new Archs4ExpressionMetadata();
	
	private String hdfExpressionFile;
	// Some data-set names we will be using.
	private static String expressionDSName = "/data/expression";
	private static String normExpressionDSName = "/data/normalized_expression";
	private static String genesDSName = "/meta/genes";
	private static String tissueDSName = "/meta/Sample_source_name_ch1";
	private static String sampleIdDSName = "/meta/Sample_geo_accession";
	private static String sampleLastUpdatedDSName = "/meta/Sample_last_update_date";
	// Eventually, we will need to *filter* the genes from the Expression file: genes not in the Correlation file will need to be excluded.

	private static Map<String,Object> expressionValuesCache = new HashMap<>();
	
	/**
	 * Contains metadata for an HDF file with ARCHS4 data.
	 * Use for the original file from ARCHS4, as well as the HDF 
	 * containing normalized expression values.
	 * @author sshorser
	 *
	 */
	class Archs4ExpressionMetadata
	{
		// Need a list of all Tissue-type names.
		private Set<String> tissueTypes = new HashSet<>();
		// A mapping of gene name to row-index in expression dataset.
		private Map<String, Integer> geneIndices = new HashMap<>();
		private Map<Integer, String> geneIndicesToNames = new HashMap<>();
		// A mapping of indices and the tissue type associated with it.
		private Map<Integer, String> indexOfTissues = new HashMap<>();
		private Map<String, List<Integer>> tissueTypeToIndex = new HashMap<>();
		// keep track of the indices for sample IDs.
		private Map<String, Integer> sampleIdToIndex = new HashMap<>();
		private Map<Integer, String> sampleIndexToID = new HashMap<>();
		private Map<String, LocalDateTime> sampleUpdateDates = new HashMap<>();
		private int numberOfSamples ;
		private int numberOfGenes ;
		public Set<String> getTissueTypes()
		{
			return tissueTypes;
		}
		public void setTissueTypes(Set<String> tissueTypes)
		{
			this.tissueTypes = tissueTypes;
		}
		public Map<String, Integer> getGeneIndices()
		{
			return geneIndices;
		}
		public void setGeneIndices(Map<String, Integer> geneIndices)
		{
			this.geneIndices = geneIndices;
		}
		public Map<Integer, String> getGeneIndicesToNames()
		{
			return geneIndicesToNames;
		}
		public void setGeneIndicesToNames(Map<Integer, String> geneIndicesToNames)
		{
			this.geneIndicesToNames = geneIndicesToNames;
		}
		public Map<Integer, String> getIndexOfTissues()
		{
			return indexOfTissues;
		}
		public void setIndexOfTissues(Map<Integer, String> indexOfTissues)
		{
			this.indexOfTissues = indexOfTissues;
		}
		public Map<String, List<Integer>> getTissueTypeToIndex()
		{
			return tissueTypeToIndex;
		}
		public void setTissueTypeToIndex(Map<String, List<Integer>> tissueTypeToIndex)
		{
			this.tissueTypeToIndex = tissueTypeToIndex;
		}
		public Map<String, Integer> getSampleIdToIndex()
		{
			return sampleIdToIndex;
		}
		public void setSampleIdToIndex(Map<String, Integer> sampleIdToIndex)
		{
			this.sampleIdToIndex = sampleIdToIndex;
		}
		public Map<Integer, String> getSampleIndexToID()
		{
			return sampleIndexToID;
		}
		public void setSampleIndexToID(Map<Integer, String> sampleIndexToID)
		{
			this.sampleIndexToID = sampleIndexToID;
		}
		public Map<String, LocalDateTime> getSampleUpdateDates()
		{
			return sampleUpdateDates;
		}
		public void setSampleUpdateDates(Map<String, LocalDateTime> sampleUpdateDates)
		{
			this.sampleUpdateDates = sampleUpdateDates;
		}
		public int getNumberOfSamples()
		{
			return numberOfSamples;
		}
		public void setNumberOfSamples(int numberOfSamples)
		{
			this.numberOfSamples = numberOfSamples;
		}
		public int getNumberOfGenes()
		{
			return numberOfGenes;
		}
		public void setNumberOfGenes(int numberOfGenes)
		{
			this.numberOfGenes = numberOfGenes;
		}
	}
	/**
	 * Used to determine what we are looking at when
	 * querying objects in an HDF file.
	 * When you call H5.H5Gget_obj_info_all on a named
	 * object in the file, the type will have a numeric value
	 * that should be in this enum.
	 * @author sshorser
	 *
	 */
	enum H5O_type
	{
		H5O_TYPE_UNKNOWN(-1), // Unknown object type
		H5O_TYPE_GROUP(0), // Object is a group
		H5O_TYPE_DATASET(1), // Object is a dataset
		H5O_TYPE_NAMED_DATATYPE(2), // Object is a named data type
		H5O_TYPE_NTYPES(3); // Number of different object types
		private static final Map<Integer, H5O_type> lookup = new HashMap<>();

		static
		{
			for (H5O_type s : EnumSet.allOf(H5O_type.class))
			{
				lookup.put(s.getCode(), s);
			}
		}

		private int code;

		H5O_type(int layout_type)
		{
			this.code = layout_type;
		}

		public int getCode()
		{
			return this.code;
		}

		public static H5O_type get(int code)
		{
			return lookup.get(code);
		}
	}
	
	Archs4ExpressionDataLoader(String fileName)
	{
		this.hdfExpressionFile = fileName;
	}

	/**
	 * Allows setting a metadata object at creation time. Useful if you are reading an expression file that does not have metadata stored in it.
	 * If you know it should match the metadata of another file, pass that metadata here.
	 * @param fileName
	 * @param md
	 */
	Archs4ExpressionDataLoader(String fileName, Archs4ExpressionMetadata md)
	{
		this.hdfExpressionFile = fileName;
		this.metadata = md;
	}
	
	Archs4ExpressionMetadata getMetadata()
	{
		return this.metadata;
	}
	
	/**
	 * Gets element coordinates for a gene (identified by index) across all samples for a tissue (identified by name).
	 * @param tissue - The name of the tissue.
	 * @param geneIndex - The index of a gene.
	 * @return An array that contains the coordinates in the dataset of the expression data for a specific gene across all samples for the specified tissue.
	 */
	private long[][] getElementCoordinatesForTissue(String tissue, int geneIndex)
	{
		List<Integer> indicesForTissue = this.metadata.getTissueTypeToIndex().get(tissue);
		long[][] elementCoords = new long[indicesForTissue.size()][2];
		
		for (int i = 0; i < indicesForTissue.size(); i++)
		{
			elementCoords[i][1] = geneIndex;
			elementCoords[i][0] = indicesForTissue.get(i);
		}
		return elementCoords;
	}
	
	/**
	 * Get expression values for all genes, for a tissue. All sample IDs for the tissue should be in a text file. 
	 * @param tissueFileName - the name of the file that contains the sample IDs for the tissue.
	 * @return a matrix of expression values. Columns are genes, rows are samples.
	 * @throws IOException
	 */
	public synchronized double[][] getExpressionValuesforTissue(Path tissueFileName) throws IOException
	{
		if (expressionValuesCache.containsKey(tissueFileName.toString()))
		{
			logger.trace("expression values found in cache for {}", tissueFileName.toString());
			return (double[][]) expressionValuesCache.get(tissueFileName.toString());
		}
		else
		{
			logger.info("Nothing in expression value cache for {}, loading it now...", tissueFileName.toString());
			List<String> sampleIds = Files.readAllLines(tissueFileName);
			List<Integer> indicesForTissue = new ArrayList<>();
			for (String sampleId : sampleIds)
			{
				indicesForTissue.add(this.metadata.getSampleIdToIndex().get(sampleId));
			}
			String tissueName = tissueFileName.getFileName().toString();
			double[][] values = getExpressionValuesBySampleIndices(indicesForTissue, tissueName);
			expressionValuesCache.put(tissueFileName.toString(), values);
			return values;
		}
	}

	/**
	 * Get expression values for all genes, for a tissue.
	 * @param tissue - The name of the tissue, as it is found in the dataset /meta/Sample_source_name_ch1
	 * @return a matrix of expression values. Columns are genes, rows are samples.
	 */
	public synchronized double[][] getExpressionValuesforTissue(String tissue)
	{
		if (expressionValuesCache.containsKey(tissue))
		{
			logger.trace("expression values found in cache for {}", tissue);
			return (double[][])expressionValuesCache.get(tissue);
		}
		else
		{
			logger.info("Nothing in expression value cache for {}, loading it now...", tissue);
			List<Integer> indicesForTissue = this.metadata.getTissueTypeToIndex().get(tissue);
			double[][] values = getExpressionValuesBySampleIndices(indicesForTissue, tissue);
			expressionValuesCache.put(tissue, values);
			return values;	
		}
	}
	
	/**
	 * Gets the indices for the samples that are associated with a tissue.
	 * @param tissueFileName
	 * @return
	 * @throws IOException
	 */
	public int[] getSampleIndicesForTissue(Path tissueFileName) throws IOException
	{
		List<String> sampleIds = Files.readAllLines(tissueFileName);
		int[] indicesForTissue = new int[sampleIds.size()];
		int i = 0;
		for (String sampleId : sampleIds)
		{
			indicesForTissue[i] = this.metadata.getSampleIdToIndex().get(sampleId);
			i ++;
		}
		return indicesForTissue;

	}
	
	/**
	 * Gets expression values based on a list of indices.
	 * @param indices
	 * @param datasubsetName
	 * @return
	 */
	public synchronized double[][] getExpressionValuesBySampleIndices(List<Integer> indices, String datasubsetName)
	{
		double[][] expressionValues = new double[this.metadata.getNumberOfGenes()][indices.size()];
		
		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		
		String datasetName = determineExpressionDSName(file_id);
		long dataset_id = H5.H5Dopen(file_id, datasetName, HDF5Constants.H5P_DEFAULT);
		long dset_space_id = H5.H5Dget_space(dataset_id);
		
		boolean isFirst = true;
		int op = HDF5Constants.H5S_SELECT_SET;
		
		List<List<Integer>> contigBlocks = new ArrayList<>();
		int maxBlockWidth = 0;
		int blockCount = 0;
		int blockWidth = 1;
		// make a list of contiguous blocks.
		for (int i = 0; i<indices.size(); i += blockWidth)
		{
			blockWidth = 1;
			int initialIndex = indices.get(i);
			int currentIndex = indices.get(i);
			// If there IS a next index, try to determine the block width.
			if (i + blockWidth < indices.size())
			{
				int nextIndex = indices.get(i + 1);
				// keep adding to the contiguous block, until we get a point where the next index is no longer a part of this block.
				while (nextIndex - currentIndex == 1 && i + blockWidth < indices.size() - 1)
				{
					blockWidth++;
					currentIndex = nextIndex;
					nextIndex = indices.get(i + blockWidth);
				}
			}
			// If there is no next index, then this block's width is 1.
			else 
			{
				blockWidth = 1;
			}
			
			maxBlockWidth = Math.max(blockWidth, maxBlockWidth);
			blockCount++;
			// A block is defined as having a starting point and a width. Create a "block" and add it to the list.
			List<Integer> currentBlock = Arrays.asList( initialIndex, blockWidth );
			contigBlocks.add(currentBlock);
		}
		logger.info("{} contiguous blocks; max block width: {}", blockCount, maxBlockWidth);
		
		// Append tissue-sample columns to the hyperslab...
		for (List<Integer> contigBlock : contigBlocks)
		{
			// Select a column (for a tissue sample) from the first gene to the last.
			long[] start = { contigBlock.get(0), 0 };
			long[] count = { 1, 1 };
			long[] block = { contigBlock.get(1), this.metadata.getNumberOfGenes() };
			int status = H5.H5Sselect_hyperslab(dset_space_id, op, start, null, count, block);
			if (status < 0)
			{
				logger.error("Selection returned an error code: {}", status);
			}
			// After the first selection, need to use H5S_SELECT_OR to combine multiple selections.
			if (isFirst)
			{
				isFirst = false;
				op = HDF5Constants.H5S_SELECT_OR;
			}
		}
		
		// Now... read the selections into a dataset.
		// DimX is how many total indices there were for a tissue.
		int dimx = indices.size();
		// DimY is how many genes there were (assume all, for all tissues).
		int dimy = this.metadata.getNumberOfGenes();

		// Now that the selections have all been made, it is FINALLY time to read the data.
		expressionValues = HDFUtils.readData(dataset_id, dset_space_id, dimx, dimy);
		H5.H5close();
		return expressionValues;
	}

	/**
	 * Gets the name of dataset containing expression data, for some file_id.
	 * If you are passing in a file_id for the original HDF file, the dataset name should be "expression".
	 * If you are passing in a file_id for the HDF file containing normalized data, the dataset name should be "normalized_expression".
	 * @param file_id
	 * @return the name of the dataset containing expression data.
	 */
	private static String determineExpressionDSName(long file_id)
	{
		int memberCount = (int) H5.H5Gn_members(file_id, "/");
		String[] oname = new String[memberCount];
		int[] otype = new int[memberCount];
		int[] ltype = new int[memberCount];
		long[] orefs = new long[memberCount];
		H5.H5Gget_obj_info_all(file_id, "/data", oname, otype, ltype, orefs, HDF5Constants.H5_INDEX_NAME);
		for (int indx = 0; indx < otype.length; indx++)
		{
			switch (H5O_type.get(otype[indx]))
			{
				case H5O_TYPE_DATASET:
					if ("expression".equals(oname[indx]))
					{
						return expressionDSName;
					}
					if ("normalized_expression".equals(oname[indx]))
					{
						return normExpressionDSName;
					}
					break;
			}
		}
		return null;
	}
	
	/**
	 * Gets a list of expression values for a gene/tissue pair, across all samples for the tissue.
	 * @param gene - A gene symbol/id.
	 * @param tissue - The name of the tissue.
	 * @return
	 */
	public double[] getExpressionValuesForGeneAndTissue(String gene, String tissue)
	{
		int geneIndex = this.metadata.getGeneIndices().get(gene);
		
		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		String dsName = determineExpressionDSName(file_id);
		long dataset_id = H5.H5Dopen(file_id, dsName, HDF5Constants.H5P_DEFAULT);
		long space_id = H5.H5Dget_space(dataset_id);

		List<Integer> indicesForTissue = this.metadata.getTissueTypeToIndex().get(tissue);

		long[][] elementCoords = getElementCoordinatesForTissue(tissue, geneIndex);
		
		int status = H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, elementCoords.length, elementCoords);
		if (status < 0)
		{
			logger.error("Error selecting elements! response code: {}",status);
			return null;
		}

		int dimx = indicesForTissue.size();
//		int dimy = 1;
		
		double[] dset_data = HDFUtils.readData(dataset_id, space_id, dimx);
		H5.H5close();

		double[] expressionValues = extractExpressionValuesFromDataset(dset_data);
		return expressionValues;
	}
	
	/**
	 * Gets expression values for a gene, if their "last updated date" is BEFORE a cutoff date.
	 * @param gene
	 * @param cutoff
	 * @return
	 */
	public double[] getDateFilteredExpressionValuesForGene(String gene, LocalDateTime cutoff)
	{
		double expressionValues[];
		double[] dset_data;
		// Get the indices for samples that are before the cutoff.
		List<String> sampleIDsToUse = this.metadata.getSampleUpdateDates().keySet().parallelStream()
													.filter(sampleID -> this.metadata.getSampleUpdateDates().get(sampleID).compareTo(cutoff) < 0)
													.collect(Collectors.toList());
		int geneIndex = this.metadata.getGeneIndices().get(gene);
		long[][] elementCoords = new long[sampleIDsToUse.size()][2];
		// Iterate over the sample IDs that passed the date filter.
		for (int i = 0; i < sampleIDsToUse.size(); i++)
		{
			int sampleIndex = this.metadata.getSampleIdToIndex().get(sampleIDsToUse.get(i));
			
			elementCoords[i][1] = geneIndex;
			elementCoords[i][0] = sampleIndex;
		}
		int dimx = sampleIDsToUse.size();
//		int dimy = 1;
		synchronized(expressionValuesCache)
		{
			dset_data = readExpressionValues(elementCoords, dimx);
		}
		expressionValues = extractExpressionValuesFromDataset(dset_data);
		return expressionValues;
	}

	/**
	 * Reads expression data from an HDF file.
	 * @param elementCoords Coordinates of the points to read.
	 * @param dimx
	 * @param dimy
	 * @return
	 */
//	private double[][] readExpressionValues(long[][] elementCoords, int dimx, int dimy)
//	{
//		double[][] dset_data;
//		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
//		String dsName = determineExpressionDSName(file_id);
//		long dataset_id = H5.H5Dopen(file_id, dsName, HDF5Constants.H5P_DEFAULT);
//		long space_id = H5.H5Dget_space(dataset_id);
//		
//		H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, elementCoords.length, elementCoords);
//		dset_data = HDFUtils.readData(dataset_id, space_id, dimx, dimy);
//		H5.H5close();
//		return dset_data;
//	}

	/**
	 * Reads expression data from an HDF file.
	 * @param elementCoords Coordinates of the points to read.
	 * @param dimx
	 * @return
	 */
	private double[] readExpressionValues(long[][] elementCoords, int dimx)
	{
		double[] dset_data;
		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		String dsName = determineExpressionDSName(file_id);
		long dataset_id = H5.H5Dopen(file_id, dsName, HDF5Constants.H5P_DEFAULT);
		long space_id = H5.H5Dget_space(dataset_id);
		
		H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, elementCoords.length, elementCoords);
		dset_data = HDFUtils.readData(dataset_id, space_id, dimx);
		H5.H5close();
		return dset_data;
	}
	
	/**
	 * Extracts expression values from a dset_data (a result from H5DRead)
	 * @param dset_data
	 * @return
	 */
	private static double[] extractExpressionValuesFromDataset(double[] dset_data)
	{
		double[] expressionValues;
		expressionValues = new double[dset_data.length];
		for (int i = 0; i < dset_data.length; i ++)
		{
//			for (int j = 0; j < dset_data[0].length; j++)
//			{
				// the input array is a 2-dimensional array because HDFUtils.readData always 
				// returns double[][], though I expect that for expression
				// data, the second dimension will always have a length of 1.
				// If there is MORE than one element for this inner loop to process, then 
				// something rather strange might be happening...
				double expressionValue = dset_data[i];
				expressionValues[i] = expressionValue;
//			}
		}
		return expressionValues;
	}
	
	/**
	 * Gets all expression values for a gene.
	 * @param gene - A gene-symbol.
	 * @return an array of expression values, as integers.
	 */
	public double[] getAllExpressionValuesForGene(String gene)
	{
		double expressionValues[];
		double[] dset_data;
		int geneIndex = this.metadata.getGeneIndices().get(gene);
		// Iterate over ALL tissues.
		long[][] elementCoords = new long[this.metadata.getIndexOfTissues().keySet().size()][2];
		for (int i = 0; i < this.metadata.getIndexOfTissues().keySet().size(); i++)
		{
			elementCoords[i][1] = geneIndex;
			elementCoords[i][0] = i;
		}
		int dimx = this.metadata.getIndexOfTissues().keySet().size();
//		int dimy = 1;
		synchronized(expressionValuesCache)
		{
			dset_data = readExpressionValues(elementCoords, dimx);
			
		}
		expressionValues = extractExpressionValuesFromDataset(dset_data);
		return expressionValues;
	}
	
	/**
	 * Loads sample indicies and IDs from the HDF file into an in-memory cache.
	 */
	public void loadSampleIndices()
	{
		this.metadata.setSampleIdToIndex(new HashMap<>());
		
		StringBuffer[] str_data = HDFUtils.readDataSet(this.hdfExpressionFile, Archs4ExpressionDataLoader.sampleIdDSName, this.metadata.getNumberOfSamples());
		StringBuffer[] str_data_dates = HDFUtils.readDataSet(this.hdfExpressionFile, Archs4ExpressionDataLoader.sampleLastUpdatedDSName, this.metadata.getNumberOfSamples());
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("MMM dd yyyy");
		
		logger.info("Number of elements: {}", str_data.length);
		for (int indx = 0; indx <  str_data.length; indx++)
		{
			String sampleId = str_data[indx].toString();
			this.metadata.getSampleIdToIndex().put(sampleId, indx);
			this.metadata.getSampleIndexToID().put(indx, sampleId);
			
			LocalDateTime ldt = LocalDate.parse( str_data_dates[indx].toString(), formatter).atStartOfDay();
			this.metadata.getSampleUpdateDates().put(sampleId, ldt);
		}
		logger.info("Number of Sample IDs loaded: {}", this.metadata.getSampleIdToIndex().size());
	}

	/**
	 * loads gene names, and their indices, from the HDF file into an in-memory cache.
	 */
	public void loadGeneNames()
	{
		this.metadata.setGeneIndices(new HashMap<>());
		StringBuffer[] str_data = HDFUtils.readDataSet(this.hdfExpressionFile, Archs4ExpressionDataLoader.genesDSName, this.metadata.getNumberOfGenes());
		logger.info("Number of elements: {}", str_data.length);
		for (int indx = 0; indx < str_data.length; indx++)
		{
			this.metadata.getGeneIndices().put(str_data[indx].toString(), indx);
			this.metadata.getGeneIndicesToNames().put(indx, str_data[indx].toString());
		}
		logger.info("Number of genes loaded: {}", this.metadata.getGeneIndices().size());
	}
	
	/**
	 * Returns a set of all tissue type names.
	 * @return
	 */
	public Set<String> getTissueTypes()
	{
		return this.metadata.getTissueTypes();
	}
	
	/**
	 * Loads tissue type names, and their indices, from the HDF file into an in-memory cache.
	 */
	public void loadTissueTypeNames()
	{
		this.metadata.setTissueTypeToIndex(new HashMap<>());
		this.metadata.setIndexOfTissues(new HashMap<>());

		StringBuffer[] str_data = HDFUtils.readDataSet(this.hdfExpressionFile,tissueDSName, this.metadata.getNumberOfSamples());
		logger.info("Number of elements: ", str_data.length);
		for (int indx = 0; indx < str_data.length; indx++)
		{
			String tissueType = str_data[indx].toString();
			this.metadata.getTissueTypes().add(tissueType);
			this.metadata.getIndexOfTissues().put(indx, tissueType);
			if (this.metadata.getTissueTypeToIndex().containsKey(tissueType))
			{
				this.metadata.getTissueTypeToIndex().get(tissueType).add(indx);
			}
			else
			{
				List<Integer> list = new ArrayList<>();
				list.add(indx);
				this.metadata.getTissueTypeToIndex().put(tissueType, list);
			}
		}

		logger.info("Number of distinct elements: {}", this.metadata.getTissueTypes().size());
	}
	
	/**
	 * Gets a mapping of tissue indices, mapped by their names.
	 * @return
	 */
	public Map<Integer, String> getIndexOfTissues()
	{
		return this.metadata.getIndexOfTissues();
	}

	/**
	 * Gets a mapping of tissue types mapped to their indices.
	 * @return
	 */
	public Map<String, List<Integer>> getTissueTypeToIndex()
	{
		return this.metadata.getTissueTypeToIndex();
	}


	/**
	 * Gets the path to the HDF file containing the gene expression data.
	 * @return
	 */
	public String getHdfExpressionFile()
	{
		return hdfExpressionFile;
	}

	/**
	 * Gets a mapping of gene indices, mapped by their names.
	 * @return
	 */
	public Map<String, Integer> getGeneIndices()
	{
		return this.metadata.getGeneIndices();
	}
	
	/**
	 * Sets the path to the HDF file containe gene expression data, and then load some data.
	 * NOTE: This method will ALSO attempt to read the file and load the number of
	 * genes and samples from the file.
	 * @param hdfExpressionFile - the full path to the HDF file.
	 */
	public void loadCounts()
	{
		long file_id = H5.H5Fopen(this.hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		String dsName = determineExpressionDSName(file_id);
		long dataset_id = H5.H5Dopen(file_id, dsName, HDF5Constants.H5P_DEFAULT);
		long space_id = H5.H5Dget_space(dataset_id);
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		this.metadata.setNumberOfSamples((int) dims[0]); // You should get ~167k here.
		this.metadata.setNumberOfGenes((int) dims[1]); // You should get ~35k here.
	}
	
	/**
	 * Loads various metadata: gene-name-to-index mapping, tissue names, tissue-index-to-names mapping,
	 * tissue-names-to-indices mapping, sample-names-to-index mapping
	 */
	public void loadMetaData()
	{
		this.loadGeneNames();
		this.loadTissueTypeNames();
		this.loadSampleIndices();
	}

	/**
	 * Gets a mapping of gene indices mapped to gene names.
	 * @return
	 */
	public Map<Integer, String> getGeneIndicesToNames()
	{
		return this.metadata.getGeneIndicesToNames();
	}

	/**
	 * Gets a mapping of sample indices mapped to sample IDs.
	 * @return
	 */
	public Map<Integer, String> getSampleIndexToID()
	{
		return this.metadata.getSampleIndexToID();
	}
}
