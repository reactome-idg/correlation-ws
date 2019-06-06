package org.reactome.idg.dao;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.List;
import java.util.Map;

import org.hibernate.FlushMode;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.exception.ConstraintViolationException;
import org.hibernate.query.NativeQuery;
import org.reactome.idg.model.GenePairCorrelation;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Isolation;
import org.springframework.transaction.annotation.Transactional;

// Maybe this should also implement org.reactome.idg.DataRepository
// OR... deprecate org.reactome.idg.DataRepository if we don't need so many different Data Repositories?

@Repository
public class GeneCorrelationDAOImpl implements GeneCorrelationDAO
{
	private Provenance currentProvenance;
	// These should be configurable...
	private int numTxOps = 0;
	private int batchSize = 100;
	// TODO: get dbname from config.
	private String dbName = "correlation_db";
	
	private Transaction addGeneTx;
	
	@Autowired
	private SessionFactory sessionFactory;

	private Session session;
	
	public GeneCorrelationDAOImpl()
	{
		// no-op constructor
	}

	public GeneCorrelationDAOImpl(Provenance p)
	{
		this.currentProvenance = p;
	}

	@Override
	@Transactional(isolation = Isolation.READ_UNCOMMITTED)
	public void loadGenePairsFromDataFile(String pathToFile)
	{
		session = sessionFactory.getCurrentSession();
		if (session == null || !session.isOpen())
		{
			this.session = sessionFactory.openSession();
		}
		// play with some tuning parameters... 
//		this.session.createSQLQuery("set global innodb_buffer_pool_size=8053063680;").executeUpdate();
//		this.session.createSQLQuery("set global innodb_io_capacity=5000;").executeUpdate();
//		this.session.createSQLQuery("set global innodb_io_capacity_max=20000;").executeUpdate();
//		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 2;").executeUpdate();
		this.session.createSQLQuery("SET global unique_checks=0;").executeUpdate();
		this.session.createSQLQuery("SET global autocommit=0;").executeUpdate();
//		this.session.createSQLQuery("SET sql_log_bin ='OFF';").executeUpdate();
		LocalDateTime start = LocalDateTime.now();
		NativeQuery<?> nq = this.session.createSQLQuery("LOAD DATA LOCAL INFILE :file INTO TABLE "+dbName+".gene_pair_correlation"
				+ " FIELDS ENCLOSED BY \"'\" LINES TERMINATED BY '\\n' "
				+ " (gene_1, gene_2, correlation_value, provenance_id);").setParameter("file", pathToFile);
		long numRows = nq.executeUpdate();
		LocalDateTime end = LocalDateTime.now();
//		this.session.createSQLQuery("SET sql_log_bin ='ON';").executeUpdate();
		this.session.createSQLQuery("SET global unique_checks=1;").executeUpdate();
//		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 1;").executeUpdate();
		System.out.println("Number of rows loaded: " + numRows + ", time duration: " + Duration.between(start, end).toString());
	}
	
	/**
	 * Adds a gene-pair correlation value to the database. NOTE: The provenance that this gene-pair must be set separately, using <code>setCurrentProvenance()</code>.
	 * @param gene1 - the first gene of the pair.
	 * @param gene2 - the second gene of the pair.
	 * @param correlationValue - the Correlation value.
	 */
	@Override
	public void addGenePair(String gene1, String gene2, double correlationValue)
	{
		// NOTE: This method is NOT @Transactional because I want to have control over when commits happen.
		// Executing a commit EVERY time a record is added is very slow. Executing a commit after some large number
		// of INSERTs speeds things up considerably.
		
			
		if (session == null || !session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		
		session.setHibernateFlushMode(FlushMode.COMMIT);

		GenePairCorrelation correlation = new GenePairCorrelation(gene1, gene2, correlationValue, this.currentProvenance);
		try
		{
			if (this.numTxOps == 0)
			{
				this.addGeneTx = session.beginTransaction();
			}
			session.save(correlation);
			this.numTxOps++;

			// commit when we've executed enough operations.
			if (this.numTxOps > this.batchSize - 1)
			{
				this.numTxOps = 0;
				addGeneTx.commit();
			}
		}
		catch (ConstraintViolationException e)
		{
			if (this.addGeneTx.isActive())
			{
				addGeneTx.rollback();
			}
			if (e.getCause().getMessage().contains("Duplicate entry") && e.getCause().getMessage().contains("idx_gene_pair_provenance"))
			{
				System.out.println("It looks like the gene-pair " + gene1 + "," + gene2 + " has already been loaded for that provenance (ID: "+this.currentProvenance.getId()+").");
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw e;
		}
		
	}

	/**
	 * Adds a Provenance object to the database.
	 */
	@Override
	@Transactional
	public Provenance addProvenance(Provenance p)
	{
		Provenance createdProvenance;
		
		session = sessionFactory.getCurrentSession();
		
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		
//			session.setHibernateFlushMode(FlushMode.COMMIT);
		// Before we try to persis this, let's make sure that it's not already there.
		@SuppressWarnings("unchecked")
		List<Provenance> results = session.createQuery("from Provenance where name = :name and url = :url and category = :cat and subcategory = :subcat")
											.setParameter("name", p.getName())
											.setParameter("url", p.getUrl())
											.setParameter("cat", p.getCategory())
											.setParameter("subcat", p.getSubcategory())
											.getResultList();
		

		if (null == results || results.size() == 0)
		{
			
			Long createdProvenanceId;
			createdProvenanceId = (Long) session.save(p);
			createdProvenance = (Provenance) session.createQuery("from Provenance where id = :id")
													.setParameter("id", createdProvenanceId)
													.getResultList().get(0);
		}
		else
		{
			System.out.println("Provenance already exists, and will not be recreated.");
			createdProvenance = results.get(0);
		}
		
		return createdProvenance;
	}
	
	/**
	 * Get all correlations for a pair of gene names.
	 */
	@Override
	public Map<Provenance, Double> queryForCorrelation(String gene1, String gene2)
	{
//		Session session = sessionFactory.getCurrentSession();
//		List<T> correlations = session.createQuery("SELECT * FROM gene_pair_correlation WHERE gene_1 = :g1 AND gene_2 = :g2 ")
//									.setParameter("g1", gene1)
//									.setParameter("g2", gene2)
//									.getResultList();
//		return correlations;
		return null;
	}

	/**
	 * Get the correlation value for a pair of genes for a given Provenance.
	 */
	@Override
	public Double queryForCorrelation(String gene1, String gene2, Provenance prov)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Provenance getCurrentProvenance()
	{
		return currentProvenance;
	}

	/**
	 * Because many gene-pairs can be inserted for a single provenance, it makes sense that "current provenance"
	 * is a part of the state of this DAO. Adding a gene pair will be done in the context of the current provenance.
	 * Use this method to set the current Provenance. 
	 * @param currentProvenance - the Provenance that all gene-pairs will be used when inserting gene-pairs.
	 */
	@Override
	public void setCurrentProvenance(Provenance currentProvenance)
	{
		this.currentProvenance = currentProvenance;
	}

	/**
	 * Batch size is the number of inserts before a COMMIT is executed.
	 * @return the current batch size for this DAO.
	 */
	public int getBatchSize()
	{
		return this.batchSize;
	}

	/**
	 * Batch size is the number of inserts before a COMMIT is executed. This function sets the batch size.
	 * @param batchSize - How many gene-pair correlation insertions should happen before issuing a COMMIT. 
	 */
	public void setBatchSize(int batchSize)
	{
		this.batchSize = batchSize;
	}

	/**
	 * How many opertaions have occurred since the last COMMIT.
	 * @return The number of operations in the current gene-pair correlation insertion transaction. 
	 */
	public int getNumTxOps()
	{
		return this.numTxOps;
	}
}
