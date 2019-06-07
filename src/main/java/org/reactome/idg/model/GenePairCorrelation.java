package org.reactome.idg.model;

import javax.persistence.Column;
import javax.persistence.ConstraintMode;
import javax.persistence.Entity;
import javax.persistence.ForeignKey;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

import org.hibernate.annotations.GenericGenerator;
import org.hibernate.annotations.Parameter;
import org.hibernate.id.enhanced.SequenceStyleGenerator;

@Entity
@Table(name = "gene_pair_correlation", indexes = { @Index(columnList = "gene_1, gene_2", unique = false),
													@Index(columnList = "gene_1,gene_2,provenance_id", unique = true, name = "idx_gene_pair_provenance")})
public class GenePairCorrelation
{
//	public GenePairCorrelation()
//	{
//
//	}
	
	public GenePairCorrelation(String gene1, String gene2, double correlationValue, Provenance provenance)
	{
		this.gene1 = gene1;
		this.gene2 = gene2;
		this.correlationValue = correlationValue;
		this.provenance = provenance;
	}
	
	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id;

	@Column(name = "gene_1", nullable = false, length = 40)
	private String gene1;
	
	@Column(name = "gene_2", nullable = false, length = 40)
	private String gene2;
	
	@Column(name = "correlation_value", nullable = false)
	private double correlationValue;
	
	@JoinColumn(name = "provenance_id", foreignKey = @ForeignKey(value = ConstraintMode.CONSTRAINT), nullable = false)
	@ManyToOne
	private Provenance provenance;

	public GenePairCorrelation(Long id)
	{
		this.id = id;
	}
 
	public Long getId()
	{
		return id;
	}
	
}