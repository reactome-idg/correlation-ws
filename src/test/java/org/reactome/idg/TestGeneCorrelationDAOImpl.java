package org.reactome.idg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.sql.SQLException;
import java.util.Map;

import javax.sql.DataSource;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.dao.GeneCorrelationDAO;
import org.reactome.idg.dao.GeneCorrelationDAOImpl;
import org.reactome.idg.dao.ProvenanceDAO;
import org.reactome.idg.model.Provenance;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

@SuppressWarnings("static-method")
public class TestGeneCorrelationDAOImpl
{
	
	@Test
	public void testLoadFileIntoDatabaseIT()
	{
		//NOTE: before running this test, drop gene_pair_correlation and provenance - do that outside of the hibernate context (hibernate will try to recreate these tables - it gets locked up if you try to drop them).
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			
			GeneCorrelationDAO dao = (GeneCorrelationDAO) context.getBean("dao");
			
			Provenance p = new Provenance();
			p.setId(1);
			p.setName("TEST");
			Provenance p1 = dao.addProvenance(p);
			assertNotNull(p1);
			dao.loadGenePairsFromDataFile("/tmp/small_data_for_idg");
		}
	}
	
	@Test
	public void testGetCorrelationIT()
	{
		// NOTE this test relies on only loading the ARCHS4 dataset.
		
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			
			GeneCorrelationDAO dao = (GeneCorrelationDAO) context.getBean("dao");
			Map<Provenance, Double> correlations = dao.getCorrelation("A1BG", "A1BG");
			assertTrue(correlations.size() > 0);
			assertTrue(correlations.keySet().size() == 1);
			
			Provenance p = correlations.keySet().stream().findFirst().orElse(null);
			
			Map<Provenance, Double> correlations1 = dao.getCorrelation("A1BG", "A1CF");
			Map<Provenance, Double> correlations2 = dao.getCorrelation("A1CF", "A1BG");
			
			Provenance p1 = correlations1.keySet().stream().findFirst().orElse(null);
			Provenance p2 = correlations2.keySet().stream().findFirst().orElse(null);
			
			// We assume that only the ACHS4 dataset has been loaded, so the provenance ID should be the same in both cases.
			// this test may need to be rewritten when multiple datasets get loaded.
			assertEquals(p1.getId(), p2.getId());
			
			assertTrue(correlations1.size() > 0);
			assertNotNull(correlations1.get(p1));

			assertTrue(correlations2.size() > 0);
			assertNotNull(correlations2.get(p2));

			System.out.println("Correlation value for A1BG,A1CF: "+correlations1.get(p1));
			System.out.println("Correlation value for A1CF,A1BG: "+correlations2.get(p2));
			assertTrue(correlations1.get(p1) == 0.311017245661844);
			
			assertEquals(correlations1.get(p1), correlations2.get(p2), 0);
			
		}
	}
	
	@Test
	public void testGetCorrelationByProvenanceIT()
	{
		// NOTE this test relies on only loading the ARCHS4 dataset.
		
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			
			GeneCorrelationDAO dao = (GeneCorrelationDAO) context.getBean("dao");
			
			ProvenanceDAO provenanceDao = (ProvenanceDAO) context.getBean("provenanceDao");
			
			Provenance provenance = provenanceDao.getProvenanceById(new Long(1));
			
			Double correlationValue = dao.getCorrelation("A1BG", "A1CF", provenance);
			
			assertEquals(correlationValue, 0.311017245661844, 0);
			
		}
	}
}
