import com.google.common.base.Joiner;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class ProteinNetwork implements Cloneable {
    private LinkedHashMap<Integer, Protein> proteinNetwork;
    public ArrayList<Double> effectorMotorRates;
    public ArrayList<Double> realisticBugParameters;
    public ArrayList<Double> tauLogSpace;
    private int lastProteinId;
    public double responseScalingFactor, bugStepSize, epsilon;
    public int howManyStepsInUnitTime;

    public ProteinNetwork() {
        proteinNetwork = new LinkedHashMap<Integer, Protein>();
        effectorMotorRates = new ArrayList<Double>();
        realisticBugParameters = new ArrayList<Double>();
        tauLogSpace = new ArrayList<Double>();

        for(int i = 1; i <= 4; i++) {
            realisticBugParameters.add(0.0);
        }
    }

    public void initialise() {
        //boolean containsKinase = false;

        for(lastProteinId = Protein.EXTRA_INTERACTIONS; lastProteinId < Parameters.initialProteinCount + Protein.EXTRA_INTERACTIONS; lastProteinId++) {
            Protein protein = new Protein(lastProteinId);

            /*if(protein.isKinase())
                containsKinase = true;*/

            proteinNetwork.put(lastProteinId, protein);
        }
        /*if(!containsKinase) {
            int proteinId = Modularity.random.nextInt(Parameters.initialProteinCount) + Protein.EXTRA_INTERACTIONS;
            proteinNetwork.get(proteinId).switchRole();
        }*/

        lastProteinId--;
        loadRates();
    }

    @SuppressWarnings("unchecked") //this is to keep java from complaining about the unchecked cast from object to hashmap in pn.proteinNetwork = (HashMap<Integer, Protein>) pn.proteinNetwork.clone();
    public Object clone() throws CloneNotSupportedException {
        ProteinNetwork pn = null;

        try {
            pn = (ProteinNetwork) super.clone();
        } catch(CloneNotSupportedException e) {
            System.out.println("clone() in ProteinNetwork: CloneNotSupportedException");
            System.exit(-47);
        }

        pn.proteinNetwork = (LinkedHashMap<Integer, Protein>) pn.proteinNetwork.clone();
        LinkedHashMap<Integer, Protein> tmp = new LinkedHashMap<Integer, Protein>();

        for(Map.Entry<Integer, Protein> entry : pn.proteinNetwork.entrySet()) {
            Protein p = (Protein) entry.getValue().clone();
            tmp.put(entry.getKey(), p);
        }
        for(Map.Entry<Integer, Protein> entry : tmp.entrySet())
            pn.proteinNetwork.put(entry.getKey(), entry.getValue());

        pn.effectorMotorRates = (ArrayList<Double>) pn.effectorMotorRates.clone();
        pn.realisticBugParameters = (ArrayList<Double>) pn.realisticBugParameters.clone();
        pn.tauLogSpace = (ArrayList<Double>) pn.tauLogSpace.clone();
        return pn;
    }

    //sensitivity
    public void setSensitivity(double sensitivity) {
        Set<Integer> receptorIds = getReceptorIds();

        for(Integer i : receptorIds)
            proteinNetwork.get(i).addInteraction(Protein.RECEPTOR_SENSITIVITY, sensitivity);
    }

    //post-mutagenesis
    public void updateInteractionsAfterDuplication(int originalProteinId, int duplicateProteinId) {
        Set<Integer> proteinsToUpdate = new HashSet<Integer>(proteinNetwork.keySet());
        proteinsToUpdate.remove(originalProteinId);
        proteinsToUpdate.remove(duplicateProteinId);

        for(Integer i : proteinsToUpdate)
            proteinNetwork.get(i).updateInteractionsAfterDuplication(originalProteinId, duplicateProteinId);
    }

    public void updateInteractionsAfterDeletion(int proteinId) {
        Set<Integer> proteinsToUpdate = new HashSet<Integer>(proteinNetwork.keySet());

        for(Integer i : proteinsToUpdate)
            proteinNetwork.get(i).updateInteractionsAfterDeletion(proteinId);
    }

    //getRandomProteinKey
    public int getRandomProteinKey(Set<Integer> keys, Set<Integer> exclude) {
        if(keys == null) {
            keys = new HashSet<Integer>(proteinNetwork.keySet()); //both keyset and hashset use pointers, they do NOT clone - but here it doesn't matter since we do not actually change the elements in the set - we optionally change the set (by excluding some elements from it), but we do not change the elements as such. BUT proteinNetwork.keySet() must NOT be used instead of HashSet<Integer>(proteinNetwork.keySet()) since excluding some elements in the former case would lead to their exclusion from the hashmap the keyset was produced from.
        }
        if(exclude != null) {
            keys.removeAll(exclude);
            if(keys.isEmpty()) return -2;
        }

        Integer[] keysArray = keys.toArray(new Integer[0]);
        return keysArray[Modularity.random.nextInt(keysArray.length)];
    }

    //getters and setters - protein network
    public int size() {
        return proteinNetwork.size();
    }

    public int increaseLastProteinId() {
        lastProteinId++;
        return lastProteinId;
    }

    public int getInteractionCount() {
        int interactionCount = 0;

        for(Protein p : proteinNetwork.values()) {
            interactionCount += p.getInteractionIds().size();

            if(p.isReceptor()) {
                interactionCount--;
            }
        }

        return interactionCount;
    }

    //getters and setters - proteins
    public Protein getProtein(int proteinId) {
        return proteinNetwork.get(proteinId);
    }

    public void addProtein(int proteinId, Protein protein) {
        proteinNetwork.put(proteinId, protein);
    }

    public void removeProtein(int proteinId) {
        proteinNetwork.remove(proteinId);
    }

    //getters and setters - receptors
    public int getReceptorId() {
        for(Protein p : proteinNetwork.values()) {
            if(p.isReceptor()) {
                return p.getProteinId();
            }
        }

        return -2; //error code for no receptors in the network
    }

    public Set<Integer> getReceptorIds() {
        Set<Integer> receptorIds = new HashSet<Integer>();

        for(Protein p : proteinNetwork.values()) {
            if(p.isReceptor())
                receptorIds.add(p.getProteinId());
        }

        return receptorIds;
    }

    public int getReceptorCount() {
        int receptorCount = 0;

        for(Protein p : proteinNetwork.values())
            if(p.isReceptor())
                receptorCount++;

        return receptorCount;
    }

    public Set<Integer> getProteinIds() {
        return proteinNetwork.keySet();
    }

    //print
    public void print() {
        String s;
        NumberFormat numberFormat = new DecimalFormat("00");

        for(Map.Entry<Integer, Protein> protein : proteinNetwork.entrySet()) {
            s = numberFormat.format(protein.getKey());
            System.out.print("Protein " + s + (protein.getValue().isKinase() ? " (kinase)     " : " (phosphatase)") + " interacts with proteins ");
            protein.getValue().print();

            if(protein.getValue().isReceptor())
                System.out.print("receptor");
            if(protein.getValue().isEffector())
                System.out.print("effector");

            System.out.println();
        }

        Joiner j = Joiner.on(", ");

        for(Double d : tauLogSpace) {
            System.out.format("%f, ", Math.exp(d));
        }

        System.out.println();
        System.out.println(j.join(effectorMotorRates) + ", " + responseScalingFactor + ", " + j.join(realisticBugParameters));
        System.out.println();
    }

    public void appendToFile(int generation) {
        NumberFormat format = new DecimalFormat("0.######E0");

        if(Parameters.whatToSimulate == 2) {
            if(Parameters.useIdealisedBug) {
                Modularity.effectorMotorRates.println(generation + "," + format.format(effectorMotorRates.get(0)) + "," + format.format(effectorMotorRates.get(1)) + "," + format.format(responseScalingFactor));
            } else if(Parameters.useRealisticBug) {
                Modularity.effectorMotorRates.println(generation + "," + format.format(effectorMotorRates.get(0)) + "," + format.format(effectorMotorRates.get(1)) + "," + format.format(realisticBugParameters.get(0)) + "," + format.format(realisticBugParameters.get(1)) + "," + format.format(realisticBugParameters.get(2)) + "," + format.format(realisticBugParameters.get(3)) + "," + format.format(epsilon));
            } else {
                Modularity.effectorMotorRates.println(generation + "," + format.format(effectorMotorRates.get(0)) + "," + format.format(effectorMotorRates.get(1)));
            }
        }

        Modularity.network.println("generation " + generation);
        Modularity.interactions.println("generation " + generation);
        Modularity.receptors.println("generation " + generation);
        Modularity.roles.println("generation " + generation);
        Modularity.effectors.println("generation " + generation);

        for(Map.Entry<Integer, Protein> protein : proteinNetwork.entrySet()) {
            Set<Integer> interactionIds = protein.getValue().getInteractionIds();
            Modularity.receptors.println(protein.getKey() + " = " + (protein.getValue().isReceptor() ? "isReceptor" : "isNotReceptor"));
            Modularity.roles.println(protein.getKey() + " = " + (protein.getValue().isKinase() ? "isKinase" : "isPhosphatase"));
            Modularity.effectors.println(protein.getKey() + " = " + (protein.getValue().isEffector() ? "isEffector" : "isNotEffector"));

            for(Integer i : interactionIds) {
                if(i == 0) {
                    Modularity.network.println(protein.getKey() + " deactivates " + protein.getKey());
                    Modularity.interactions.println(protein.getKey() + " (deactivates) " + protein.getKey() + " = " + format.format(protein.getValue().getInteraction(i)));
                }
                else if(i == 1) {
                    Modularity.network.println(protein.getKey() + " activates " + protein.getKey());
                    Modularity.interactions.println(protein.getKey() + " (activates) " + protein.getKey() + " = " + format.format(protein.getValue().getInteraction(i)));
                }
                else if(i == 2) {
                    Modularity.network.println("F activates " + protein.getKey());
                    Modularity.interactions.println("F (activates) " + protein.getKey() + " = " + format.format(protein.getValue().getInteraction(i)));
                }
                else {
                    Modularity.network.println(protein.getKey() + (protein.getValue().isKinase() ? " activates " : " deactivates ") + i);
                    Modularity.interactions.println(protein.getKey() + (protein.getValue().isKinase() ? " (activates) " : " (deactivates) ") + i + " = " + format.format(protein.getValue().getInteraction(i)));
                }
            }
        }

        Modularity.network.println();
        Modularity.interactions.println();
        Modularity.receptors.println();
        Modularity.roles.println();
        Modularity.effectors.println();
    }

    private double roundToTwoDecimalPlaces(double d) {
        NumberFormat twoDecimalPlaces = new DecimalFormat("#.##");
        return Double.valueOf(twoDecimalPlaces.format(d));
    }

    public void loadRates() {
        effectorMotorRates.add(Parameters.effectorMotorAssociationRate);
        effectorMotorRates.add(Parameters.effectorMotorDissociationRate);
    }

    public void loadResponseScalingFactor() {
        responseScalingFactor = Parameters.responseScalingFactor;
    }

    public void loadRealisticBug() {
        realisticBugParameters.set(0, Parameters.realisticBugConstant);
        realisticBugParameters.set(1, Parameters.realisticBugLinear);
        realisticBugParameters.set(2, Parameters.realisticBugQuadratic);
        realisticBugParameters.set(3, Parameters.realisticBugMemoryLength);
        epsilon = Parameters.epsilon;

        if(Parameters.mutationRateOriginal[14] > 0) {
            calculateABFromEpsilonSensitivity();
        }

        if(Parameters.calculateLinearFromQuadratic) {
            calculateLinear();
        }

        if(Parameters.calculateBugStepSize) {
            calculateBugStepSize();
        } else {
            calculateBugStepSize(Parameters.idealisedBugStepSize);

            if(bugStepSize != Parameters.idealisedBugStepSize) {
                System.out.println("bugStepSize is different from Parameters.idealisedBugStepSize! bugStepSize = " + bugStepSize + ", Parameters.idealisedBugStepSize = " + Parameters.idealisedBugStepSize);
            }
        }

        if(Modularity.whatTodo == 3 && Parameters.whatToSimulate == 1 && Parameters.stimulusCourse[5]) {
            Parameters.stimulusCourseDuration = realisticBugParameters.get(3) * Parameters.memoryLengthFactor + 1;

            if(Parameters.stimulusCourseDuration > 3000) {
                Parameters.stimulusCourseDuration = 3000;
            }
        }
    }

    public void calculateLinear() {
        if(Parameters.useConstantAndLinearForAdaptive) {
            realisticBugParameters.set(0, -(realisticBugParameters.get(1) / (1 + Parameters.epsilon)) * (1 - Parameters.epsilon));
        } else {
            realisticBugParameters.set(1, -2 * realisticBugParameters.get(2));
        }
    }

    public void calculateBugStepSize(double d) {
        howManyStepsInUnitTime = (int) Math.ceil(Parameters.stepSize / d);

        if(howManyStepsInUnitTime % 2 == 1) {
            howManyStepsInUnitTime++;
        }

        bugStepSize = Parameters.stepSize / howManyStepsInUnitTime;
    }

    public void calculateBugStepSize() {
        calculateBugStepSize(Math.min(Math.min(1 / Parameters.stochasticRate, Parameters.worldLength / (Parameters.numberOfModes * Bug.VELOCITY_MAX)), realisticBugParameters.get(3)) / 20);
    }

    public void fillTauLogSpace() {
        tauLogSpace.add(Math.log(effectorMotorRates.get(0)));
        tauLogSpace.add(Math.log(effectorMotorRates.get(1)));
        tauLogSpace.add(Math.log(responseScalingFactor));

        if(!Parameters.useNewWayOfMutatingAAndB) {
            tauLogSpace.add(Math.log(realisticBugParameters.get(0)));
            tauLogSpace.add(Math.log(realisticBugParameters.get(1)));
        } else {
            tauLogSpace.add(Math.signum(realisticBugParameters.get(0)) * (Math.log(Math.abs(realisticBugParameters.get(0))) - Parameters.minimumPossibleDecimalValue));
            tauLogSpace.add(Math.signum(realisticBugParameters.get(1)) * (Math.log(Math.abs(realisticBugParameters.get(1))) - Parameters.minimumPossibleDecimalValue));
        }

        tauLogSpace.add(Math.log(realisticBugParameters.get(2)));
        tauLogSpace.add(Math.log(realisticBugParameters.get(3)));
    }

    public void calculateABFromEpsilonSensitivity() {
        if(epsilon <= 0) {
            realisticBugParameters.set(0, Parameters.sensitivity);
            realisticBugParameters.set(1, -realisticBugParameters.get(0) * (1 + epsilon) / (1 - epsilon));
        } else {
            realisticBugParameters.set(1, -Parameters.sensitivity);
            realisticBugParameters.set(0, -realisticBugParameters.get(1) * (1 - epsilon) / (1 + epsilon));
        }
    }

    public void setLastProteinId() {
        for(Integer i : proteinNetwork.keySet()) {
            if(lastProteinId < i) lastProteinId = i;
        }
    }
}
