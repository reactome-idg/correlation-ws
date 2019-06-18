package org.reactome.idg.loader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
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

public class H5ExpressionDataLoader
{
	private static final Logger logger = LogManager.getLogger();
	// Need a list of all Tissue-type names.
	private static Set<String> tissueTypes = new HashSet<>();
	// A mapping of gene name to row-index in expression dataset.
	private static Map<String, Integer> geneIndices = new HashMap<>();
	// A mapping of indices and the tissue type associated with it.
	private static Map<Integer, String> indexOfTissues = new HashMap<>();
	private static Map<String, List<Integer>> tissueTypeToIndex = new HashMap<>();

	// TODO: NUM_ROWS can probably be queried directly from the file, no need to hard-code it.
	private static final int NUM_SAMPLES=167726;
	private static String FILENAME = "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5";
	String[] Xname = new String[1];
	// Some data-set names we will be using.
	private final static String expressionDSName = "/data/expression";
	private final static String genesDSName = "/meta/genes";
	private final static String tissueDSName = "/meta/Sample_source_name_ch1";
	// TODO: NUM_GENES_IN_FILE can probably be queried from the file. No need to hard code it.
	private static final int NUM_GENES_IN_FILE = 35238;
	
	
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
	

	public static int[][] getExpressionValuesforTissue(String tissue)
	{
		List<Integer> indicesForTissue = tissueTypeToIndex.get(tissue);
		logger.info("number of samples for tissue ({}): {}", tissue, indicesForTissue.size());
		logger.info("Tissue indices: {}", indicesForTissue.toString());
		int[][] expressionValues = new int[NUM_GENES_IN_FILE][indicesForTissue.size()];
		
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
		long dset_space_id = H5.H5Dget_space(dataset_id);
		
		// I think I can iterate over all ranges, and ADD them to a selection and then do one single select at the end. To try that later today...
		// The only problem is that the output doesn't seem right... :/
		// Append selections for each gene
		int geneCount = 0;
		boolean isFirst = true;
		int op = HDF5Constants.H5S_SELECT_SET;
		
		List<List<Integer>> contigBlocks = new ArrayList<>();
		int maxBlockWidth = 0;
		int blockCount = 0;
		int blockWidth = 1;
		// make a list of contiguous blocks.
		for (int i = 0; i<indicesForTissue.size(); i += blockWidth)
		{
			blockWidth = 1;
			int initialIndex = indicesForTissue.get(i);
			int currentIndex = indicesForTissue.get(i);
			if (i + blockWidth < indicesForTissue.size())
			{
				int nextIndex = indicesForTissue.get(i + 1);
				// keep adding to the contiguous block, until we get a point where the next index is no longer a part of this block.
				while (nextIndex - currentIndex == 1 && i + blockWidth < indicesForTissue.size() - 1)
				{
					blockWidth++;
					currentIndex = nextIndex;
					nextIndex = indicesForTissue.get(i + blockWidth);
				}
			}
			else 
			{
				blockWidth = 1;
			}
			
			maxBlockWidth = Math.max(blockWidth, maxBlockWidth);
			blockCount++;
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
			long[] block = { contigBlock.get(1), NUM_GENES_IN_FILE };
			int status = H5.H5Sselect_hyperslab(dset_space_id, op, start, null, count, block);
			if (status < 0)
			{
				logger.error("Selection returned an error code: {}", status);
			}
			if (isFirst)
			{
				isFirst = false;
				op = HDF5Constants.H5S_SELECT_OR;
			}
		}
		
		// Now... read the selections into a dataset.
		// DimX is how many total indices there were for a tissue.
		int dimx = indicesForTissue.size();
		// DimY is how many genes there were (assume all, for all tissues).
		int dimy = NUM_GENES_IN_FILE;		

		expressionValues = readData(dataset_id, dset_space_id, dimx, dimy);
		H5.H5close();
		return expressionValues;
	}
	

	/**
	 * Gets a list of expression values for a gene/tissue pair.
	 * @param gene
	 * @param tissue
	 * @return
	 */
	public static int[] getExpressionValuesForGeneAndTissue(String gene, String tissue)
	{
		int geneIndex = geneIndices.get(gene);
		
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);

		List<Integer> indicesForTissue = tissueTypeToIndex.get(tissue);

		long space_id = H5.H5Dget_space(dataset_id);

		long[][] elementCoords = getElementCoordinatesForTissue(tissue, geneIndex);
		
		int status = H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, elementCoords.length, elementCoords);
		if (status < 0)
		{
			logger.error("Error selecting elements! response code: {}",status);
			return null;
		}
		long[] start = new long[2];
		long[] end = new long[2];
		H5.H5Sget_select_bounds(space_id, start, end);
		int dimx = indicesForTissue.size();
		int dimy = 1;
		
		int[][] dset_data = readData(dataset_id, space_id, dimx, dimy);
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


	private static int[][] readData(long dataset_id, long space_id, int dimx, int dimy)
	{
		int[][] dset_data = new int[dimx][dimy];
		
		long[] dims = new long[2];
		dims[0] = dimx;
		dims[1] = dimy;
		long memspace_id = H5.H5Screate_simple(2, dims, null);
		long type_id = H5.H5Dget_type(dataset_id);
		H5.H5Dread(dataset_id, type_id, memspace_id, space_id, HDF5Constants.H5P_DEFAULT, dset_data);
		H5.H5Sclose(memspace_id);
		H5.H5Tclose(type_id);

		return dset_data;
	}
	
	public static void loadGeneNames()
	{
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, genesDSName, HDF5Constants.H5P_DEFAULT);
		long type_id = H5.H5Dget_type(dataset_id);
		long dataWidth = H5.H5Tget_size(type_id);
		long space_id = H5.H5Dget_space(dataset_id);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		byte[][] dset_data = new byte[NUM_GENES_IN_FILE][(int) dataWidth];
		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
		byte[] tempbuf = new byte[(int) dataWidth];
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			for (int jndx = 0; jndx < dataWidth; jndx++)
			{
				tempbuf[jndx] = dset_data[indx][jndx];
			}
			str_data[indx] = new StringBuffer(new String(tempbuf).trim());
		}
		H5.H5Dclose(dataset_id);
		H5.H5Tclose(type_id);
		H5.H5Sclose(space_id);
		H5.H5Fclose(file_id);
		logger.info("Number of elements: {}",maxdims[0]);
		for (int indx = 0; indx < maxdims[0]; indx++)
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
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, tissueDSName, HDF5Constants.H5P_DEFAULT);
		long type_id = H5.H5Dget_type(dataset_id);
		long dataWidth = H5.H5Tget_size(type_id);
		long space_id = H5.H5Dget_space(dataset_id);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		byte[][] dset_data = new byte[NUM_SAMPLES][(int) dataWidth];
		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
		byte[] tempbuf = new byte[(int) dataWidth];
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			for (int jndx = 0; jndx < dataWidth; jndx++)
			{
				tempbuf[jndx] = dset_data[indx][jndx];
			}
			str_data[indx] = new StringBuffer(new String(tempbuf).trim());
		}
		H5.H5Dclose(dataset_id);
		H5.H5Tclose(type_id);
		H5.H5Sclose(space_id);
		H5.H5Fclose(file_id);
		logger.info("Number of elements: ", maxdims[0]);
		for (int indx = 0; indx < maxdims[0]; indx++)
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
		Map<Integer,Integer> histogram = new HashMap<>();
		for (String tisssue : tissueTypeToIndex.keySet())
		{
			int size = tissueTypeToIndex.get(tisssue).size();
			if (histogram.containsKey(size))
			{
				histogram.put(size, histogram.get(size)+1 );
			}
			else
			{
				histogram.put(size, 1);
			}
		}
		for (Integer size : histogram.keySet().stream().sorted().collect(Collectors.toList()))
		{
			logger.info("Samples / tissue: {} ; # different tissues: {}", size, histogram.get(size));
		}
		
		logger.info("Number of distinct elements: {}", tissueTypes.size());
	}
	
	public static void getExpressionValuesForGene(String geneId)
	{
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		
		if (file_id >= 0)
		{
			long dataset_id = H5.H5Dopen(file_id, genesDSName, HDF5Constants.H5P_DEFAULT);
			
			long type_id = H5.H5Dget_type(dataset_id);
			System.out.println(type_id);
			long dataWidth = H5.H5Tget_size(type_id);
			long space_id = H5.H5Dget_space(dataset_id);
			long[] dims = new long[2];
			long[] maxdims = new long[2];
			H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
			System.out.println(dims[0]);
			System.out.println(maxdims[0]);
			byte[][] dset_data = new byte[NUM_SAMPLES][(int) dataWidth];
			System.out.println(dataWidth);
			if (dataset_id >= 0)
			{
				long dcpl_id = H5.H5Dget_create_plist(dataset_id);
				if (dcpl_id >= 0)
				{
					if (dataset_id >= 0)
					{
						StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
						long memtype_id = H5.H5Tcopy(HDF5Constants.H5T_C_S1);
						H5.H5Tset_size(memtype_id, dataWidth);
						H5.H5Dread(dataset_id, memtype_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
						
						byte[] tempbuf = new byte[(int) dataWidth];
						for (int indx = 0; indx < maxdims[0]; indx++)
						{
							for (int jndx = 0; jndx < dataWidth; jndx++)
							{
								tempbuf[jndx] = dset_data[indx][jndx];
							}
							str_data[indx] = new StringBuffer(new String(tempbuf).trim());
						}
						
						for (int indx = 0; indx < maxdims[0]; indx++)
						{
							logger.debug("{} [{}]: {}", genesDSName, indx, str_data[indx]);
						}
					}
				}
			}
			
			dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
			type_id = H5.H5Dget_type(dataset_id);
//			System.out.println(type_id);
//			System.out.println(H5.H5Tget_tag(type_id));
		}
	}

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


	public static Map<Integer, String> getIndexOfTissues() {
		return indexOfTissues;
	}


	public static Map<String, List<Integer>> getTissueTypeToIndex() {
		return tissueTypeToIndex;
	}
	
}
