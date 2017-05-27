import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Protein implements Cloneable {
    private HashMap<Integer, Double> proteinInteractions;
    private int proteinId;
    public static final int EXTRA_INTERACTIONS = 3; //rate of passive deactivation, rate of passive activation, sensitivity of receptor to stimulus
    public static final int RECEPTOR_SENSITIVITY = 2;
    private boolean isReceptor = false, isEffector = false, isKinase;

    public Protein(int proteinId) {
        this.proteinId = proteinId;
        proteinInteractions = new HashMap<Integer, Double>();

        if((Modularity.whatTodo == 1 || Modularity.whatTodo == 2) && Parameters.loadNetwork == -1) {
            isKinase = Parameters.proteinRoles[proteinId - EXTRA_INTERACTIONS];
        } else if(Modularity.whatTodo == 3 || Modularity.whatTodo == 4 || Modularity.whatTodo == 9 || Parameters.loadNetwork != -1) {
            isKinase = true;
        }

        if(proteinId == EXTRA_INTERACTIONS) {
            isReceptor = true;
        }
        if(proteinId == EXTRA_INTERACTIONS + 1) {
            isEffector = true;
        }

        //isKinase = Modularity.random.nextBoolean();
    }

    public Protein(int proteinId, boolean isKinase) {
        this.proteinId = proteinId;
        proteinInteractions = new HashMap<Integer, Double>();
        this.isKinase = isKinase;
    }

    @SuppressWarnings("unchecked") //this is to keep java from complaining about the unchecked cast from object to hashmap in p.proteinInteractions = (HashMap<Integer, Double>) p.proteinInteractions.clone();
    public Object clone() throws CloneNotSupportedException {
        Protein p = null;

        try {
            p = (Protein) super.clone();
        } catch(CloneNotSupportedException e) {
            System.out.println("clone() in Protein: CloneNotSupportedException");
            System.exit(-47);
        }

        p.proteinInteractions = (HashMap<Integer, Double>) p.proteinInteractions.clone();
        return p;
    }

    //getters
    public boolean isReceptor() {
        return isReceptor;
    }

    public boolean isEffector() {
        return isEffector;
    }

    public boolean isKinase() {
        return isKinase;
    }

    public int getProteinId() {
        return proteinId;
    }

    public boolean containsInteraction(int proteinId) {
        return proteinInteractions.containsKey(proteinId);
    }

    public double getInteraction(int proteinId) {
        return proteinInteractions.get(proteinId);
    }

    public Set<Integer> getInteractionIds() {
        return proteinInteractions.keySet();
    }

    //setters
    public void setReceptor(boolean b) {
        isReceptor = b;
    }

    public void setEffector(boolean b) {
        isEffector = b;
    }

    //mutagenesis - mutate interaction
    public void addInteraction(int proteinId, double interactionStrength) {
        proteinInteractions.put(proteinId, interactionStrength);
    }

    //mutagenesis - duplication
    public void setProteinId(int proteinId) {
        this.proteinId = proteinId;
    }

    public void updateInteractionsAfterDuplication(int originalProteinId, int duplicateProteinId) {
        if(proteinInteractions.containsKey(originalProteinId))
            proteinInteractions.put(duplicateProteinId, proteinInteractions.get(originalProteinId));
    }

    //mutagenesis - deletion
    public void deleteInteraction(int proteinId) {
        proteinInteractions.remove(proteinId);
    }

    public void updateInteractionsAfterDeletion(int proteinId) {
        if(proteinInteractions.containsKey(proteinId))
            proteinInteractions.remove(proteinId);
    }

    //mutagenesis - switch kinase/phosphatase
    public void switchRole() {
        isKinase = !isKinase;
    }

    //print
    public void print() {
        String s;
        NumberFormat interactionId = new DecimalFormat("00");
        NumberFormat interactionStrength = new DecimalFormat("0.000");

        for(Map.Entry<Integer, Double> interaction : proteinInteractions.entrySet()) {
            s = interactionId.format(interaction.getKey());
            System.out.print(s + " (");
            s = interactionStrength.format(interaction.getValue());
            System.out.print(s + "), ");
        }
    }
}
