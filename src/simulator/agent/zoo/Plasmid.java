package simulator.agent.zoo;

import java.math.BigInteger;
import java.util.ArrayList;

import idyno.SimTimer;
import simulator.agent.InfoAgent;
import simulator.agent.LocatedAgent;
import utils.ExtraMath;
import utils.LogFile;

/**
 * \brief TODO
 * 
 * <p>NOTE: This uses simulation time, not agentGrid time (there is no
 * such thing in iDynoMiCS v1.x), so I (Rob) strongly recommend keeping 
 * agentTimeStep equal to global time step when using this class.</p> 
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 *
 */
public class Plasmid extends InfoAgent
{
	/**
	 * Number of copies of this plasmid in the current host.
	 */
	protected int _copyNumber;
	
	/**
	 * Simulation time (h) at which this plasmid last donated to another host.
	 */
	protected double _tLastDonated;
	
	/**
	 * Simulation time (h) at which this plasmid was received by its host.
	 */
	protected double _tReceived;
	
	/**
	 * Tally variable for the number of potential recipients this plasmid can
	 * scan in an agentTimeStep.
	 */
	protected double _testTally;
	
	/**
	 * Rate (in units 1/h) which which this plasmid's pilus can currently scan
	 * neighbouring cells to find a potential recipient.
	 */
	protected double _scanRate;
	/**
	 * Qian 10.2016:
	 * Number of horizontal transfers (conjugation) 
	 */
	protected int _numHT;
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public Plasmid()
	{
		super();
		_speciesParam = new PlasmidParam();
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		Plasmid out = (Plasmid) super.clone();
		out._speciesParam = this._speciesParam;
		out._generation = this._generation;
		String tempString = new String(this._genealogy.toString());
		out._genealogy = new StringBuffer(tempString);
		out._testTally = this._testTally;
		out._copyNumber = this._copyNumber;
		out._tLastDonated = this._tLastDonated;
		out._tReceived = this._tReceived;
		out._numHT = this._numHT;
		out._scanRate = this._scanRate;
		return out;
	}
	
	@Override
	public Plasmid sendNewAgent() throws CloneNotSupportedException
	{
		Plasmid baby = (Plasmid) this.clone();
		baby.reset();
		return baby;
	}
	
	/**
	 * Clone this plasmid and register the clone (in the species population).
	 */
	@Override
	public void createNewAgent()
	{
		try
		{
			Plasmid baby = this.sendNewAgent();
			baby.registerBirth();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "Plasmid.createNewAgent()");
		}
	}

	
	
	/*************************************************************************
	 * BASIC METHODS
	 ************************************************************************/
	
	/**
	 * \brief Determines if the given name is the same as this Plasmid's 
	 * species name.
	 * 
	 * @param name String to check (case-sensitive)
	 * @return boolean: true if name is species name, false if different.
	 */
	public boolean isSpeciesName(String name)
	{
		return name.equals(this._species.speciesName);
	}
	
	/**
	 * \brief Set all of this Plasmid's details to default values. 
	 */
	public void reset()
	{
		this._generation = 0;
		this._genealogy = new StringBuffer("0");
		this._copyNumber = getSpeciesParam().copyNumDefault;
		this._tLastDonated = -Double.MAX_VALUE;
		this._tReceived = -Double.MAX_VALUE;
		this._testTally = 0.0;
		this._numHT = 0;
		this._scanRate = 0.0;
	}
	
	@Override
	public PlasmidParam getSpeciesParam()
	{
		return (PlasmidParam) this._speciesParam;
	}
	
	/**
	 * No need to register birth on the agent grid - only need to tell the
	 * species.
	 */
	@Override
	public void registerBirth()
	{
		this._species.notifyBirth();
	}
	
	/**
	 * No need to register death on the agent grid - only need to tell the
	 * species.
	 */
	public void die()
	{
		this._species.notifyDeath();
	}
	
	public double getPilusRange()
	{
		return this.getSpeciesParam().pilusLength;
	}
	
	public int getPlasmidID()
	{
		return this.getSpeciesParam().plasmidID;
	}
	
	public int getCopyNumber()
	{
		return this._copyNumber;
	}
	
	public double getTimeRecieved()
	{
		return this._tReceived;
	}
	
	public double getTimeLastDonated()
	{
		return this._tLastDonated;
	}
	
	public int getNumHT()
	{
		return this._numHT;
	}
	
	/**
	 * \brief Set the copy number, time last received, and time last donated 
	 * for this Plasmid.
	 * 
	 * <p>Should be used during PlamsidBac initialisation from results
	 * file.</p>
	 * 
	 * @param copyNum int copy number.
	 * @param tReceived double time last received.
	 * @param tLastDonated double time last donated.
	 */
	public void setDetails(int copyNum, double tReceived, double tLastDonated)
	{
		this._copyNumber = copyNum;
		this._tReceived = tReceived;
		this._tLastDonated = tLastDonated;
	}
	
	/**
	 * \brief Provides a list of all the reactions encoded by this Plasmid.
	 * 
	 * @return ArrayList<Integer> of reaction indices.
	 */
	public ArrayList<Integer> getReactionsEncoded()
	{
		return this.getSpeciesParam().reactionsEncoded;
	}
	
	/*************************************************************************
	 * HIGHER METHODS
	 ************************************************************************/
	
	public boolean isHostCompatible(PlasmidBac targetRecipient)
	{
		return getSpeciesParam().isHostCompatible(targetRecipient);	    
	}

	/**
	 * \brief Given a potential recipient, determine whether this agent could
	 * receive a copy of this plasmid.
	 * 
	 * @param targetRecipient Any PlasmidBac.
	 * @return Whether or not this PlasmidBac could host this plasmid.
	 */
	public boolean isCompatible(PlasmidBac targetRecipient)
	{
		return getSpeciesParam().isCompatible(targetRecipient);
	}
	
	/**
	 * \brief Whether this Plasmid is ready to attempt donation yet.
	 * 
	 * <p>It may not conjugate if it has been lost, or if it donated or was
	 * received too recently.</p>
	 * 
	 * @return boolean: true if conjugation possible, false if not.
	 */
	public Boolean isReadyToConjugate()
	{
		/*
		 * First check if the plasmid has sufficient copies.
		 */
		if ( this._copyNumber < 1 )
			return false;
		/*
		 * Now check timings: cannot conjugate if given/received a plasmid too
		 * recently.
		 */
		double now = SimTimer.getCurrentTime();
		if ( now < this._tLastDonated + getSpeciesParam().donationLag )
			return false;
		return ( now >= this._tReceived + getSpeciesParam().receptionLag );
	}
	
	/**
	 * \brief Update this plasmid's scan rate based on the growth tone of its
	 * host.
	 * 
	 * <p>See PlasmidBac.getScaledTone() and .growthTone() for more info.</p>
	 * 
	 * @param scaledTone Scaled growth tone of the PlasmidBac hosting this 
	 * plasmid.
	 */
	public void updateScanRate(double scaledTone)
	{
		this._scanRate = this.getSpeciesParam().scanSpeed * scaledTone;
		this._testTally +=  this._scanRate * _agentGrid.AGENTTIMESTEP;
		LogFile.writeLog("scaledTone " + scaledTone + " testTally " + _testTally);
	}
	
	/**
	 * \brief Try to send a copy of this Plasmid to a target LocatedAgent. Used both in chemostat and biofilm simulations.
	 * 
	 * @param aTarget LocatedAgent that may receive a copy of this Plasmid.
	 */
	public void tryToSendPlasmid(LocatedAgent aTarget)
	{
		// static reference for counting number of tries
		PlasmidBac._numTry++;
//		LogFile.writeLog("try to send Plasmid");
		/*
		 * We're looking at a target, so decrement the _testTally to reflect this.
		 */
		this._testTally -= 1.0; // also decremented in chemostat simulations but without consequences as _testTally not evaluated
		/*
		 * Unless this is a PlasmidBac, there can be no donation.
		 */
		if ( ! ( aTarget instanceof PlasmidBac) )
			return;
		PlasmidBac aPB = (PlasmidBac) aTarget;
		/*
		 * If the target is incompatible, or this plasmid fails a proficiency
		 * test, there is no donation.
		 */
		if ( ! ( this.isCompatible(aPB) &&  this.testProficiency() ) )
			return;	
		/*
		 * Donation is successful, so make the new plasmid, give it to the
		 * target, and update the old plasmid.
		 */
		try
		{
			Plasmid baby = this.sendNewAgent();
			baby.registerBirth();
			aPB.welcomePlasmid(baby);
			baby._copyNumber = this._copyNumber;
			this._numHT++; //Qian 10.2016: update the number of horizontal transfers for this plasmid (donor).
			PlasmidBac._numTrans++;
//			LogFile.writeLog("numHT " + this.getBirthday() + " " + _numHT);
			this._tLastDonated = baby._tReceived = SimTimer.getCurrentTime();
			baby._testTally = 0.0; // jan: was this._testTally = 0.0; but a donor should be able to conjugate several times per timestep if _scanRate * AGENTTIMESTEP is > 1
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "Plasmid.tryToSendPlasmid()");
		}
	}
	
	/**
	 * \brief Check whether this plasmid still has time to scan in this
	 * timestep. canScan() is not called in a chemostat simulation, so _testTally is effectively ignored
	 * 
	 * @return boolean: true if there is time for another scan, false if not.
	 */
	public boolean canScan()
	{
		return ( this._testTally >= 1.0 );
	}
	
	/**
	 * \brief Only transfer the plasmid with a certain probability.
	 * transferProficiency is a probability between 0 and 1 and specified in the protocol file. 
	 * 
	 * @return boolean: true if donation may proceed, false if not.
	 */
	public boolean testProficiency()
	{
		return ExtraMath.getUniRandDbl() < 
										getSpeciesParam().transferProficiency;
	}
	
	/**
	 * \brief This Plasmid should have been recently created through cell
	 * division of its host: check to see if it has been lost in the process.
	 */
	public void applySegregation()
	{
		if ( ExtraMath.getUniRandDbl() < getSpeciesParam().lossProbability )
			this._copyNumber = 0;
	}
	
}