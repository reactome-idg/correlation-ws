package org.reactome.idg.loader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;

public class Archs4ExpressionDataLoader
{
	private static final Logger logger = LogManager.getLogger();
	// Need a list of all Tissue-type names.
	private static Set<String> tissueTypes = new HashSet<>();
	// A mapping of gene name to row-index in expression dataset.
	private static Map<String, Integer> geneIndices = new HashMap<>();
	// A mapping of indices and the tissue type associated with it.
	private static Map<Integer, String> indexOfTissues = new HashMap<>();
	private static Map<String, List<Integer>> tissueTypeToIndex = new HashMap<>();
	// keep track of the indices for sample IDs.
	private static Map<String, Integer> sampleIdToIndex = new HashMap<>();

	private static int numberOfSamples ;
	private static String hdfExpressionFile;
	// Some data-set names we will be using.
	private final static String expressionDSName = "/data/expression";
	private final static String genesDSName = "/meta/genes";
	private final static String tissueDSName = "/meta/Sample_source_name_ch1";
	private final static String sampleIdDSName = "/meta/Sample_geo_accession";
	// Eventually, we will need to *filter* the genes from the Expression file: genes not in the Correlation file will need to be excluded.
	private static int numberOfGenes ;

	private Archs4ExpressionDataLoader()
	{
		// private constructor to prevent instantiation. The HDF library does NOT allow concurrent access to HDF files 
		// so having multiple instances of this class could cause problems if they all try to access the file simultaneously.
		// So everything here must be static!
	}

	/**
	 * Gets element coordinates for a gene (identified by index) across all samples for a tissue (identified by name).
	 * @param tissue - The name of the tissue.
	 * @param geneIndex - The index of a gene.
	 * @return An array that contains the coordinates in the dataset of the expression data for a specific gene across all samples for the specified tissue.
	 */
	private static long[][] getElementCoordinatesForTissue(String tissue, int geneIndex)
	{
		List<Integer> indicesForTissue = tissueTypeToIndex.get(tissue);
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
	public static int[][] getExpressionValuesforTissue(Path tissueFileName) throws IOException
	{
		List<String> sampleIds = Files.readAllLines(tissueFileName);
		List<Integer> indicesForTissue = new ArrayList<>();
		for (String sampleId : sampleIds)
		{
			indicesForTissue.add(sampleIdToIndex.get(sampleId));
		}
		String tissueName = tissueFileName.getFileName().toString();
		return getExpressionValuesByIndices(indicesForTissue, tissueName);
		
	}

	/**
	 * Get expression values for all genes, for a tissue.
	 * @param tissue - The name of the tissue, as it is found in the dataset /meta/Sample_source_name_ch1
	 * @return a matrix of expression values. Columns are genes, rows are samples.
	 */
	public static int[][] getExpressionValuesforTissue(String tissue)
	{
		List<Integer> indicesForTissue = tissueTypeToIndex.get(tissue);
		return getExpressionValuesByIndices(indicesForTissue, tissue);
	}
	
	/**
	 * Gets expression values based on a list of indices.
	 * @param indices
	 * @param datasubsetName
	 * @return
	 */
	private static int[][] getExpressionValuesByIndices(List<Integer> indices, String datasubsetName)
	{
//		List<Integer> indicesForTissue = sampleIdToIndex.get(tissue);
		logger.info("number of samples for tissue ({}): {}", datasubsetName, indices.size());
		logger.info("Tissue indices: {}", indices.toString());
		int[][] expressionValues = new int[numberOfGenes][indices.size()];
		
		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
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
			long[] block = { contigBlock.get(1), numberOfGenes };
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
		int dimy = numberOfGenes;		

		// Now that the selections have all been made, it is FINALLY time to read the data.
		expressionValues = HDFUtils.readData(dataset_id, dset_space_id, dimx, dimy);
		H5.H5close();
		return expressionValues;
	}
	
	/**
	 * Gets a list of expression values for a gene/tissue pair, across all samples for the tissue.
	 * @param gene - A gene symbol/id.
	 * @param tissue - The name of the tissue.
	 * @return
	 */
	public static int[] getExpressionValuesForGeneAndTissue(String gene, String tissue)
	{
		int geneIndex = geneIndices.get(gene);
		
		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
		long space_id = H5.H5Dget_space(dataset_id);

		List<Integer> indicesForTissue = tissueTypeToIndex.get(tissue);

		long[][] elementCoords = getElementCoordinatesForTissue(tissue, geneIndex);
		
		int status = H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, elementCoords.length, elementCoords);
		if (status < 0)
		{
			logger.error("Error selecting elements! response code: {}",status);
			return null;
		}
//		long[] start = new long[2];
//		long[] end = new long[2];
//		H5.H5Sget_select_bounds(space_id, start, end);
		int dimx = indicesForTissue.size();
		int dimy = 1;
		
		int[][] dset_data = HDFUtils.readData(dataset_id, space_id, dimx, dimy);
		int[] expressionValues = new int[dset_data.length];
		for (int i = 0; i < dset_data.length; i ++)
		{
			for (int j = 0; j < dset_data[0].length; j++)
			{
				int expressionValue = dset_data[i][j];
				expressionValues[i] = expressionValue;
			}
		}
		H5.H5close();
		return expressionValues;
	}
	
	public static void loadSampleIndices()
	{
		Archs4ExpressionDataLoader.sampleIdToIndex = new HashMap<>();
		
		StringBuffer[] str_data = HDFUtils.readDataSet(Archs4ExpressionDataLoader.hdfExpressionFile, Archs4ExpressionDataLoader.sampleIdDSName, Archs4ExpressionDataLoader.numberOfSamples);
		logger.info("Number of elements: {}", str_data.length);
		for (int indx = 0; indx <  str_data.length; indx++)
		{
			Archs4ExpressionDataLoader.sampleIdToIndex.put(str_data[indx].toString(), indx);
		}
		logger.info("Number of Sample IDs loaded: {}", Archs4ExpressionDataLoader.sampleIdToIndex.size());
	}

	public static void loadGeneNames()
	{
		Archs4ExpressionDataLoader.geneIndices = new HashMap<>();
		StringBuffer[] str_data = HDFUtils.readDataSet(Archs4ExpressionDataLoader.hdfExpressionFile, Archs4ExpressionDataLoader.genesDSName, Archs4ExpressionDataLoader.numberOfGenes);
		logger.info("Number of elements: {}", str_data.length);
		for (int indx = 0; indx < str_data.length; indx++)
		{
			geneIndices.put(str_data[indx].toString(), indx);
		}
		logger.info("Number of genes loaded: {}", geneIndices.size());
	}
	
	public static Set<String> getTissueTypes()
	{
		return tissueTypes;
	}
	
	public static void loadTissueTypeNames()
	{
		Archs4ExpressionDataLoader.tissueTypeToIndex = new HashMap<>();
		Archs4ExpressionDataLoader.indexOfTissues = new HashMap<>();

		StringBuffer[] str_data = HDFUtils.readDataSet(Archs4ExpressionDataLoader.hdfExpressionFile,tissueDSName, numberOfSamples);
		logger.info("Number of elements: ", str_data.length);
		for (int indx = 0; indx < str_data.length; indx++)
		{
			String tissueType = str_data[indx].toString();
			tissueTypes.add(tissueType);
			indexOfTissues.put(indx, tissueType);
			if (tissueTypeToIndex.containsKey(tissueType))
			{
				tissueTypeToIndex.get(tissueType).add(indx);
			}
			else
			{
				List<Integer> list = new ArrayList<>();
				list.add(indx);
				tissueTypeToIndex.put(tissueType, list);
			}
		}
		
//		tissueTypeToIndex.keySet().parallelStream().forEach(t -> {if (tissueTypeToIndex.get(t).size() <= 5) {logger.info(t);} }) ;
//		System.out.println(
//		tissueTypeToIndex.keySet().parallelStream().map( t -> tissueTypeToIndex.get(t).size()).max(Integer::compareTo)
//		);
//		tissueTypeToIndex.keySet().parallelStream().forEach(t -> {if (tissueTypeToIndex.get(t).size() == 8231) {logger.info(t);} }) ;


		// size is the key, the *number* with that size is the value.
//		Map<Integer,Integer> histogram = new HashMap<>();
//		for (String tisssue : tissueTypeToIndex.keySet())
//		{
//			int size = tissueTypeToIndex.get(tisssue).size();
//			if (histogram.containsKey(size))
//			{
//				histogram.put(size, histogram.get(size)+1 );
//			}
//			else
//			{
//				histogram.put(size, 1);
//			}
//		}
//		for (Integer size : histogram.keySet().stream().sorted().collect(Collectors.toList()))
//		{
//			logger.info("Samples / tissue: {} ; # different tissues: {}", size, histogram.get(size));
//		}
		
		logger.info("Number of distinct elements: {}", tissueTypes.size());
	}
	
//	public static void getExpressionValuesForGene(String geneId)
//	{
//		long file_id = H5.H5Fopen(hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
//		
//		if (file_id >= 0)
//		{
//			long dataset_id = H5.H5Dopen(file_id, genesDSName, HDF5Constants.H5P_DEFAULT);
//			
//			long type_id = H5.H5Dget_type(dataset_id);
//			long dataWidth = H5.H5Tget_size(type_id);
//			long space_id = H5.H5Dget_space(dataset_id);
//			long[] dims = new long[2];
//			long[] maxdims = new long[2];
//			H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
//			byte[][] dset_data = new byte[numberOfSamples][(int) dataWidth];
//			if (dataset_id >= 0)
//			{
//				long dcpl_id = H5.H5Dget_create_plist(dataset_id);
//				if (dcpl_id >= 0)
//				{
//					if (dataset_id >= 0)
//					{
//						StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
//						long memtype_id = H5.H5Tcopy(HDF5Constants.H5T_C_S1);
//						H5.H5Tset_size(memtype_id, dataWidth);
//						H5.H5Dread(dataset_id, memtype_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
//						
//						byte[] tempbuf = new byte[(int) dataWidth];
//						for (int indx = 0; indx < maxdims[0]; indx++)
//						{
//							for (int jndx = 0; jndx < dataWidth; jndx++)
//							{
//								tempbuf[jndx] = dset_data[indx][jndx];
//							}
//							str_data[indx] = new StringBuffer(new String(tempbuf).trim());
//						}
//						
//						for (int indx = 0; indx < maxdims[0]; indx++)
//						{
//							logger.info("{} [{}]: {}", genesDSName, indx, str_data[indx]);
//						}
//					}
//				}
//			}
//			
//			dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
//			type_id = H5.H5Dget_type(dataset_id);
//		}
//	}

//	public static Map<String, String> getTissueToSamples()
//	{
//		Map<String, String> tissueSampleMapping = new HashMap<>();
//		
//		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
//		long dataset_id = H5.H5Dopen(file_id, "/meta/Sample_geo_accession", HDF5Constants.H5P_DEFAULT);
//		long type_id = H5.H5Dget_type(dataset_id);
//		long dataWidth = H5.H5Tget_size(type_id);
//		long space_id = H5.H5Dget_space(dataset_id);
//		long[] dims = new long[2];
//		long[] maxdims = new long[2];
//		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
//		byte[][] dset_data = new byte[NUM_SAMPLES][(int) dataWidth];
//		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
//		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
//		byte[] tempbuf = new byte[(int) dataWidth];
//		for (int indx = 0; indx < maxdims[0]; indx++)
//		{
//			for (int jndx = 0; jndx < dataWidth; jndx++)
//			{
//				tempbuf[jndx] = dset_data[indx][jndx];
//			}
//			str_data[indx] = new StringBuffer(new String(tempbuf).trim());
//		}
//		H5.H5Dclose(dataset_id);
//		H5.H5Tclose(type_id);
//		H5.H5Sclose(space_id);
//		H5.H5Fclose(file_id);
//		logger.info("Number of elements: ", maxdims[0]);
//		for (int indx = 0; indx < maxdims[0]; indx++)
//		{
//			String sampleId = str_data[indx].toString();
//			tissueSampleMapping.put(indexOfTissues.get(indx), sampleId);
//		}
//		
//		return tissueSampleMapping;
//	}


	public static Map<Integer, String> getIndexOfTissues()
	{
		return indexOfTissues;
	}


	public static Map<String, List<Integer>> getTissueTypeToIndex()
	{
		return tissueTypeToIndex;
	}


	public static String getHdfExpressionFile()
	{
		return hdfExpressionFile;
	}

	public static Map<String, Integer> getGeneIndices()
	{
		return geneIndices;
	}
	
	/**
	 * Sets the path to the HDF file containe gene expression data, and then load some data.
	 * NOTE: This method will ALSO attempt to read the file and load the number of
	 * genes and samples from the file.
	 * @param hdfExpressionFile - the full path to the HDF file.
	 */
	public static void setHdfExpressionFileAndLoadCounts(String hdfExpressionFile)
	{
		Archs4ExpressionDataLoader.hdfExpressionFile = hdfExpressionFile;
		long file_id = H5.H5Fopen(Archs4ExpressionDataLoader.hdfExpressionFile, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
		long space_id = H5.H5Dget_space(dataset_id);
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		numberOfSamples = (int) dims[0]; // You should get ~167k here.
		numberOfGenes = (int) dims[1]; // You should get ~35k here.
	}
	
	/**
	 * Loads various metadata: gene-name-to-index mapping, tissue names, tissue-index-to-names mapping,
	 * tissue-names-to-indices mapping, sample-names-to-index mapping
	 */
	public static void loadMetaData()
	{
		Archs4ExpressionDataLoader.loadGeneNames();
		Archs4ExpressionDataLoader.loadTissueTypeNames();
		Archs4ExpressionDataLoader.loadSampleIndices();
	}
}
