import java.util.Set;

public class Dynamics implements DerivnFunction {
    private ProteinNetwork proteinNetwork;
    private double input;

    public Dynamics(ProteinNetwork proteinNetwork, double input) {
        this.proteinNetwork = proteinNetwork;
        this.input = input;
    }

    public double[] derivn(double x, double[] activityLevels) {
		double[] slope = new double[activityLevels.length];
        Set<Integer> proteinIds = proteinNetwork.getProteinIds();

        for(Integer i : proteinIds)
            slope[i - Protein.EXTRA_INTERACTIONS] = 0.0; //initialise slope

        for(Integer i : proteinIds) { //take a protein
            for(Integer j : proteinNetwork.getProtein(i).getInteractionIds()) { //and go through its interactions
                if(j == Protein.RECEPTOR_SENSITIVITY) { //sensitivity
                    slope[i - Protein.EXTRA_INTERACTIONS] = slope[i - Protein.EXTRA_INTERACTIONS] + proteinNetwork.getProtein(i).getInteraction(j) * input * (1 - activityLevels[i - Protein.EXTRA_INTERACTIONS]);
                }
                else if(j == 0) { //passive deactivation
                    slope[i - Protein.EXTRA_INTERACTIONS] = slope[i - Protein.EXTRA_INTERACTIONS] - proteinNetwork.getProtein(i).getInteraction(j) * activityLevels[i - Protein.EXTRA_INTERACTIONS];
                }
                else if(j == 1) { //passive activation
                    slope[i - Protein.EXTRA_INTERACTIONS] = slope[i - Protein.EXTRA_INTERACTIONS] + proteinNetwork.getProtein(i).getInteraction(j) * (1 - activityLevels[i - Protein.EXTRA_INTERACTIONS]);
                }
                else { //interactions with other proteins
                    if(proteinNetwork.getProtein(i).isKinase()) {
                        slope[j - Protein.EXTRA_INTERACTIONS] = slope[j - Protein.EXTRA_INTERACTIONS] + proteinNetwork.getProtein(i).getInteraction(j) * activityLevels[i - Protein.EXTRA_INTERACTIONS] * (1 - activityLevels[j - Protein.EXTRA_INTERACTIONS]);
                    }
                    else {
                        slope[j - Protein.EXTRA_INTERACTIONS] = slope[j - Protein.EXTRA_INTERACTIONS] - proteinNetwork.getProtein(i).getInteraction(j) * activityLevels[i - Protein.EXTRA_INTERACTIONS] * activityLevels[j - Protein.EXTRA_INTERACTIONS];
                    }
                }
            }
        }

		return slope;
	}
}
