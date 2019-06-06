package org.reactome.idg;

import static org.junit.Assert.assertNotNull;

import java.sql.SQLException;

import javax.sql.DataSource;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.dao.GeneCorrelationDAO;
import org.reactome.idg.dao.GeneCorrelationDAOImpl;
import org.reactome.idg.model.Provenance;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class TestGeneCorrelationDAOImpl
{
	@SuppressWarnings("static-method")
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
}
