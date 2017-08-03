package simulator.agent.zoo;

import java.util.LinkedList;
import idyno.SimTimer;

/**
 * PlasmidMemory keeps track of the times a particular type of plasmid enters a particular host and gets lost again, potentially entering again,
 * so it is a list of times of entry and departure for one particular host.
 * Plasmids are distinguished by their plasmidID
 * This is used for calculating fitness cost evolution.
 * 
 * @author Akvile Zemgulyte (akvilezemgulyte@gmail.com)
 *
 * */

public class PlasmidMemory extends Object {
	// the same kind of plasmid can enter, get lost, and enter again etc so we need a list of times
	private LinkedList<Double> tPlasmidReceived;
	private LinkedList<Double> tPlasmidLost;
	private Integer plasmidID;
	
	public PlasmidMemory(Integer plasmidID){
		
		this.plasmidID = plasmidID;
		this.tPlasmidReceived=new LinkedList<>();
		this.tPlasmidLost=new LinkedList<>();
	}
	
	public void addTPlasmidReceived (double timeReceived){
		tPlasmidReceived.addLast(timeReceived);
	}
	
	public void addTPlasmidLost (double timeLost){
		tPlasmidLost.addLast(timeLost);
	}
		
	public LinkedList<Double> timeSpentInHost(){
		LinkedList<Double>timeInHost=new LinkedList<>();
		for (int i=0;i<tPlasmidLost.size(); i++)
		{
			timeInHost.addLast(tPlasmidLost.get(i)-tPlasmidReceived.get(i)); 
		}
		if (tPlasmidLost.size()<tPlasmidReceived.size()){
			timeInHost.addLast(SimTimer.getCurrentTime()-tPlasmidReceived.getLast());
		}
		
		return timeInHost;
	}
	
	public LinkedList<Double> timeSpentLost(){
		LinkedList<Double>timeBeingLost=new LinkedList<>();
		for (int i=1;i<tPlasmidLost.size(); i++)
		{
			timeBeingLost.addLast(tPlasmidReceived.get(i)-tPlasmidLost.get(i-1)); 
		}
		return timeBeingLost;
	}
	
	public Integer getPlasmidID(){
		return plasmidID;
	}
	public double calcTimeSpentInSameHost(){
		double timeInSameHost=timeSpentInHost().getLast();
		return timeInSameHost;
	}
	
}