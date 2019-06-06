package org.reactome.idg.model;

import java.io.Serializable;
import java.util.concurrent.atomic.AtomicLong;

import org.hibernate.HibernateException;
import org.hibernate.engine.spi.SessionImplementor;
import org.hibernate.engine.spi.SharedSessionContractImplementor;
import org.hibernate.id.enhanced.SequenceStyleGenerator;

public class AssignedSequenceStyleGenerator extends SequenceStyleGenerator
{
	private static AtomicLong seqValue = new AtomicLong(0);
	
	@Override
	public Serializable generate(SharedSessionContractImplementor session, Object obj) throws HibernateException 
	{
		if (obj instanceof Identifiable)
		{
			Identifiable<?> identifiable = (Identifiable<?>) obj;
			Serializable id = identifiable.getId();
			if (id != null)
			{
				return id;
			}
			else
			{
				synchronized(seqValue)
				{
					Long l = seqValue.incrementAndGet();
//					System.out.println("sequence generated: "+l.toString());
					return l;
				}
			}
		}
		return super.generate(session, obj);
	}
}
