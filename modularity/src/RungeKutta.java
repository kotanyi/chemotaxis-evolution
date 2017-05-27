import java.util.*;

public class RungeKutta implements Cloneable {
    private double stepSizeThisStep, stepSizeNextStep;
    private HashMap<Integer, Double> activityLevels;

    public RungeKutta(Set<Integer> proteinIds) {
        stepSizeNextStep = Parameters.advanceTimeByDefault; //this is needed for cashKarp
        activityLevels = new HashMap<Integer, Double>();

        if((Modularity.whatTodo == 1 || Modularity.whatTodo == 2 || Modularity.whatTodo == 4) && Parameters.loadNetwork == -1) {
            for(int i = Protein.EXTRA_INTERACTIONS; i < Protein.EXTRA_INTERACTIONS + Parameters.initialProteinCount; i++)
                activityLevels.put(i, 0.0);
        } else if(Modularity.whatTodo == 3 || Modularity.whatTodo == 9 || Parameters.loadNetwork != -1) {
            for(Integer i : proteinIds) {
                activityLevels.put(i, 0.0);
            }
        } else if(Modularity.whatTodo == 5) {
            activityLevels.put(0, 1.0);
        }
    }

    @SuppressWarnings("unchecked") //this is to keep java from complaining about the unchecked cast from object to hashmap in (HashMap<Integer, Double>) rk.activityLevels.clone();
    public Object clone() throws CloneNotSupportedException {
        RungeKutta rk = null;

        try {
            rk = (RungeKutta) super.clone();
        } catch(CloneNotSupportedException e) {
            System.out.println("clone() in RungeKutta: CloneNotSupportedException");
            System.exit(-47);
        }

        rk.activityLevels = (HashMap<Integer, Double>) rk.activityLevels.clone();
        return rk;
    }

    public void fourthOrder(ProteinNetwork proteinNetwork, double input) {
        HashMap<Integer, Double> k1 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k2 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k3 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k4 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> tmp = new HashMap<Integer, Double>();
        HashMap<Integer, Double> slope = new HashMap<Integer, Double>();

        slope = calculateSlope(proteinNetwork, activityLevels, input, slope);
        for(Map.Entry<Integer, Double> entry : slope.entrySet())
            k1.put(entry.getKey(), Parameters.advanceTimeByDefault * entry.getValue());

        for(Map.Entry<Integer, Double> entry : k1.entrySet())
            tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue()/2);
        slope = calculateSlope(proteinNetwork, tmp, input, slope);
        for(Map.Entry<Integer, Double> entry : slope.entrySet())
            k2.put(entry.getKey(), Parameters.advanceTimeByDefault * entry.getValue());

        for(Map.Entry<Integer, Double> entry : k2.entrySet())
            tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue()/2);
        slope = calculateSlope(proteinNetwork, tmp, input, slope);
        for(Map.Entry<Integer, Double> entry : slope.entrySet())
            k3.put(entry.getKey(), Parameters.advanceTimeByDefault * entry.getValue());

        for(Map.Entry<Integer, Double> entry : k3.entrySet())
            tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue());
        slope = calculateSlope(proteinNetwork, tmp, input, slope);
        for(Map.Entry<Integer, Double> entry : slope.entrySet())
            k4.put(entry.getKey(), Parameters.advanceTimeByDefault * entry.getValue());

        for(Map.Entry<Integer, Double> entry : k1.entrySet())
            tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue()/6 + k2.get(entry.getKey())/3 + k3.get(entry.getKey())/3 + k4.get(entry.getKey())/6);

        stepSizeThisStep = stepSizeNextStep = Parameters.advanceTimeByDefault; //one might be tempted to remove this and put it in the constructor. BUT, bear in mind that when runSimulation is run, cashKarp will be used most of the time. if, then, fourthOrder will be used at some point, stepSizeThisStep and stepSizeNextStep will NOT be changed by it (= fourthOrder) and will remain the same as set by the last run of cashKarp.
        activityLevels = tmp;
    }

    public boolean cashKarp(ProteinNetwork proteinNetwork, double input) {
        double safety = 0.9, error, maxError, delta, tolerance, stepSize = stepSizeNextStep;
        HashMap<Integer, Double> k1 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k2 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k3 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k4 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k5 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k6 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> tmp = new HashMap<Integer, Double>();
        HashMap<Integer, Double> fourthOrderAccurate = new HashMap<Integer, Double>();
        HashMap<Integer, Double> fifthOrderAccurate = new HashMap<Integer, Double>();
        HashMap<Integer, Double> slope = new HashMap<Integer, Double>();

        for(int i = 0; i <= Parameters.maxNumberOfIterations; i++) {
            //System.out.println("Adaptive RK: iteration #" + i);
            slope = calculateSlope(proteinNetwork, activityLevels, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k1.put(entry.getKey(), stepSize * entry.getValue());

            for(Map.Entry<Integer, Double> entry : k1.entrySet())
                tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue() / 5);
            slope = calculateSlope(proteinNetwork, tmp, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k2.put(entry.getKey(), stepSize * entry.getValue());

            for(Map.Entry<Integer, Double> entry : k1.entrySet())
                tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + (entry.getValue()*3 + k2.get(entry.getKey())*9) / 40);
            slope = calculateSlope(proteinNetwork, tmp, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k3.put(entry.getKey(), stepSize * entry.getValue());

            for(Map.Entry<Integer, Double> entry : k1.entrySet())
                tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + (entry.getValue()*3 - k2.get(entry.getKey())*9 + k3.get(entry.getKey())*12) / 10);
            slope = calculateSlope(proteinNetwork, tmp, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k4.put(entry.getKey(), stepSize * entry.getValue());

            for(Map.Entry<Integer, Double> entry : k1.entrySet())
                tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) - entry.getValue()*11/54 + k2.get(entry.getKey())*5/2 - k3.get(entry.getKey())*70/27 + k4.get(entry.getKey())*35/27);
            slope = calculateSlope(proteinNetwork, tmp, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k5.put(entry.getKey(), stepSize * entry.getValue());

            for(Map.Entry<Integer, Double> entry : k1.entrySet())
                tmp.put(entry.getKey(), activityLevels.get(entry.getKey()) + entry.getValue()*1631/55296 + k2.get(entry.getKey())*175/512 + k3.get(entry.getKey())*575/13824 + k4.get(entry.getKey())*44275/110592 + k5.get(entry.getKey())*253/4096);
            slope = calculateSlope(proteinNetwork, tmp, input, slope);
            for(Map.Entry<Integer, Double> entry : slope.entrySet())
                k6.put(entry.getKey(), stepSize * entry.getValue());

            maxError = 0;

            for(Map.Entry<Integer, Double> entry : activityLevels.entrySet()) {
                fourthOrderAccurate.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey()) * 2825 / 27648 + k3.get(entry.getKey()) * 18575 / 48384 + k4.get(entry.getKey()) * 13525 / 55296 + k5.get(entry.getKey()) * 277 / 14336 + k6.get(entry.getKey()) / 4);
                fifthOrderAccurate.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey()) * 37 / 378 + k3.get(entry.getKey()) * 250 / 621 + k4.get(entry.getKey()) * 125 / 594 + k6.get(entry.getKey()) * 512 / 1771);
                error = Math.abs(fifthOrderAccurate.get(entry.getKey()) - fourthOrderAccurate.get(entry.getKey()));
                tolerance = Math.abs(fourthOrderAccurate.get(entry.getKey())) * Parameters.relativeTolerance + Parameters.absoluteTolerance;

                if(tolerance == 0) {
                    maxError = 2; //if tolerance == 0, error/tolerance = positive infinity and is therefore surely greater than the current maxError (see how maxError should be calculated in the else below), so we can just arbitrarily set it to 2 (see what maxError is later used for).
                }
                else {
                    maxError = Math.max(maxError, error/tolerance);
                }

                if(fourthOrderAccurate.get(entry.getKey()) > 1 || fourthOrderAccurate.get(entry.getKey()) < 0) {
                    maxError = 2;
                }
            }

            if(maxError <= 1) {
                for(Integer j : fourthOrderAccurate.keySet()) {
                    if(fourthOrderAccurate.get(j) < 0 || fourthOrderAccurate.get(j) > 1 || fifthOrderAccurate.get(j) < 0 || fifthOrderAccurate.get(j) > 1) {
                        System.out.println("Adaptive RK: protein activity is higher/lower than 1/0.");
                    }
                }

                stepSizeThisStep = stepSize;
                delta = safety * Math.pow(maxError, -0.2);

                if(delta > 4) {
                    stepSize *= 4;
                }
                else if(delta > 1) {
                    stepSize *= delta;
                }

                stepSizeNextStep = stepSize;
                activityLevels = fourthOrderAccurate;

                return true;
            }
            else {
                delta = safety * Math.pow(maxError, -0.25);

                if(delta < 0.1) {
                    stepSize *= 0.1;
                }
                else {
                    stepSize *= delta;
                }
            }
        }

        return false;
    }

    public boolean cashKarp(ProteinNetwork proteinNetwork, World world, Bug bug, boolean isEquilibrating) {
        double safety = 0.9, error, maxError, delta, tolerance, stepSize = stepSizeNextStep, timeNextStep, baseline, time = world.getTime();
        HashMap<Integer, Double> k1 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k2 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k3 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k4 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k5 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> k6 = new HashMap<Integer, Double>();
        HashMap<Integer, Double> tmp = new HashMap<Integer, Double>();
        HashMap<Integer, Double> fourthOrderAccurate = new HashMap<Integer, Double>();
        HashMap<Integer, Double> fifthOrderAccurate = new HashMap<Integer, Double>();
        HashMap<Integer, Double> slope = new HashMap<Integer, Double>();
        List<Double> inputs = new ArrayList<Double>();
        Set<Map.Entry<Integer, Double>> entrySet = activityLevels.entrySet();

        if(Parameters.whatToSimulate == 2 && isEquilibrating) {
            if(Parameters.useStochasticFood) {
                baseline = Parameters.stochasticBaseline;
            } else {
                baseline = Parameters.gaussianBaselineConcentration;
            }

            for(int j = 0; j < 6; j++) {
                inputs.add(baseline);
            }
        }

        for(int i = 0; i <= Parameters.maxNumberOfIterations; i++) {
            //System.out.println("Adaptive RK: iteration #" + i);
            if(Parameters.whatToSimulate == 2 && !isEquilibrating) {
                if(Parameters.useStochasticFood) {
                    inputs = getStochasticFood(bug.getPosition(), bug.isTumbling(), bug.getVelocity(), bug.getDriftVelocity(), stepSize, world);
                } else {
                    inputs = World.getInputGaussians(bug.getPosition(), bug.isTumbling(), bug.getVelocity(), bug.getDriftVelocity(), stepSize);
                }
            }

            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, activityLevels, World.getInput(time), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, activityLevels, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, activityLevels, inputs.get(0), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k1.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            for(Map.Entry<Integer, Double> entry : entrySet)
                tmp.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey()) / 5);
            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, World.getInput(time + stepSize / 5), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, tmp, inputs.get(1), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k2.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            for(Map.Entry<Integer, Double> entry : entrySet)
                tmp.put(entry.getKey(), entry.getValue() + (k1.get(entry.getKey())*3 + k2.get(entry.getKey())*9) / 40);
            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, World.getInput(time + stepSize * 3 / 10), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, tmp, inputs.get(2), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k3.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            for(Map.Entry<Integer, Double> entry : entrySet)
                tmp.put(entry.getKey(), entry.getValue() + (k1.get(entry.getKey())*3 - k2.get(entry.getKey())*9 + k3.get(entry.getKey())*12) / 10);
            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, World.getInput(time + stepSize * 3 / 5), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, tmp, inputs.get(3), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k4.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            for(Map.Entry<Integer, Double> entry : entrySet)
                tmp.put(entry.getKey(), entry.getValue() - k1.get(entry.getKey())*11/54 + k2.get(entry.getKey())*5/2 - k3.get(entry.getKey())*70/27 + k4.get(entry.getKey())*35/27);
            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, World.getInput(time + stepSize), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, tmp, inputs.get(4), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k5.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            for(Map.Entry<Integer, Double> entry : entrySet)
                tmp.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey())*1631/55296 + k2.get(entry.getKey())*175/512 + k3.get(entry.getKey())*575/13824 + k4.get(entry.getKey())*44275/110592 + k5.get(entry.getKey())*253/4096);
            if(Parameters.whatToSimulate == 1 && !isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, World.getInput(time + stepSize * 7 / 8), slope); }
            if(Parameters.whatToSimulate == 1 && isEquilibrating) { slope = calculateSlope(proteinNetwork, tmp, 0, slope); }
            if(Parameters.whatToSimulate == 2) { slope = calculateSlope(proteinNetwork, tmp, inputs.get(5), slope); }
            for(Map.Entry<Integer, Double> entry : entrySet)
                k6.put(entry.getKey(), stepSize * slope.get(entry.getKey()));

            maxError = 0;

            for(Map.Entry<Integer, Double> entry : entrySet) {
                fourthOrderAccurate.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey()) * 2825 / 27648 + k3.get(entry.getKey()) * 18575 / 48384 + k4.get(entry.getKey()) * 13525 / 55296 + k5.get(entry.getKey()) * 277 / 14336 + k6.get(entry.getKey()) / 4);
                fifthOrderAccurate.put(entry.getKey(), entry.getValue() + k1.get(entry.getKey()) * 37 / 378 + k3.get(entry.getKey()) * 250 / 621 + k4.get(entry.getKey()) * 125 / 594 + k6.get(entry.getKey()) * 512 / 1771);
                error = Math.abs(fifthOrderAccurate.get(entry.getKey()) - fourthOrderAccurate.get(entry.getKey()));
                tolerance = Math.abs(fourthOrderAccurate.get(entry.getKey())) * Parameters.relativeTolerance + Parameters.absoluteTolerance;

                if(tolerance == 0) {
                    maxError = 2; //if tolerance == 0, error/tolerance = positive infinity and is therefore surely greater than the current maxError (see how maxError should be calculated in the else below), so we can just arbitrarily set it to 2 (see what maxError is later used for).
                }
                else {
                    maxError = Math.max(maxError, error/tolerance);
                }

                if(fourthOrderAccurate.get(entry.getKey()) > 1 || fourthOrderAccurate.get(entry.getKey()) < 0 || fifthOrderAccurate.get(entry.getKey()) > 1 || fifthOrderAccurate.get(entry.getKey()) < 0) {
                    maxError = 2;
                }
            }

            if(maxError <= 1) {
                for(Integer j : fourthOrderAccurate.keySet()) {
                    if(fourthOrderAccurate.get(j) < 0 || fourthOrderAccurate.get(j) > 1 || fifthOrderAccurate.get(j) < 0 || fifthOrderAccurate.get(j) > 1) {
                        System.out.println("Adaptive RK: protein activity is higher/lower than 1/0.");
                    }
                }

                stepSizeThisStep = stepSize;
                delta = safety * Math.pow(maxError, -0.2);

                if(delta > 4) {
                    stepSize *= 4;
                }
                else if(delta > 1) {
                    stepSize *= delta;
                }

                timeNextStep = Parameters.stepSize * (Math.floor((time + stepSizeThisStep) / Parameters.stepSize) + 1);

                if(time + stepSizeThisStep + stepSize > timeNextStep) {
                    stepSize = timeNextStep - (time + stepSizeThisStep);
                }

                stepSizeNextStep = stepSize;
                activityLevels = fourthOrderAccurate;

                return true;
            }
            else {
                delta = safety * Math.pow(maxError, -0.25);

                if(delta < 0.1) {
                    stepSize *= 0.1;
                }
                else {
                    stepSize *= delta;
                }
            }
        }

        return false;
    }

    private List<Double> getStochasticFood(double position, boolean isTumbling, double velocity, double driftVelocity, double stepSize, World world) {
        List<Double> food = new ArrayList<Double>();
        double time = world.getTime();

        if(isTumbling) {
            velocity = driftVelocity;
        } else {
            velocity += driftVelocity;
        }

        food.add(world.stochasticFoodDistribution(position, World.getIndexOfNearestStochasticTime(time)));
        food.add(world.stochasticFoodDistribution(position + velocity * stepSize / 5, World.getIndexOfNearestStochasticTime(time + stepSize / 5)));
        food.add(world.stochasticFoodDistribution(position + velocity * stepSize * 3 / 10, World.getIndexOfNearestStochasticTime(time + stepSize * 3 / 10)));
        food.add(world.stochasticFoodDistribution(position + velocity * stepSize * 3 / 5, World.getIndexOfNearestStochasticTime(time + stepSize * 3 / 5)));
        food.add(world.stochasticFoodDistribution(position + velocity * stepSize, World.getIndexOfNearestStochasticTime(time + stepSize)));
        food.add(world.stochasticFoodDistribution(position + velocity * stepSize * 7 / 8, World.getIndexOfNearestStochasticTime(time + stepSize * 7 / 8)));
        return food;
    }

    private HashMap<Integer, Double> calculateSlope(ProteinNetwork proteinNetwork, HashMap<Integer, Double> activityLevels, double input, HashMap<Integer, Double> slope) {
        Set<Integer> proteinIds = proteinNetwork.getProteinIds();
        double interaction, activityLevelsValue;
        Protein protein;

        for(Integer i : proteinIds)
            slope.put(i, 0.0); //initialise slope

        for(Integer i : proteinIds) { //take a protein
            activityLevelsValue = activityLevels.get(i);
            protein = proteinNetwork.getProtein(i);

            for(Integer j : protein.getInteractionIds()) { //and go through its interactions
                interaction = protein.getInteraction(j);
                
                if(j == Protein.RECEPTOR_SENSITIVITY) { //sensitivity
                    slope.put(i, slope.get(i) + interaction * input * (1 - activityLevelsValue));
                }
                else if(j == 0) { //passive deactivation
                    slope.put(i, slope.get(i) - interaction * activityLevelsValue);
                }
                else if(j == 1) { //passive activation
                    slope.put(i, slope.get(i) + interaction * (1 - activityLevelsValue));
                }
                else { //interactions with other proteins
                    if(proteinNetwork.getProtein(i).isKinase()) {
                        slope.put(j, slope.get(j) + interaction * activityLevelsValue * (1 - activityLevels.get(j)));
                    }
                    else {
                        slope.put(j, slope.get(j) - interaction * activityLevelsValue * activityLevels.get(j));
                    }
                }
            }
        }

		return slope;
	}

    //update activityLevels after duplication/deletion
    public void addProtein(int proteinId) {
        activityLevels.put(proteinId, 0.0);
    }

    public void deleteProtein(int proteinId) {
        activityLevels.remove(proteinId);
    }

    //getters & setters
    public void resetStepSizeNextStep() {
        stepSizeNextStep = Parameters.advanceTimeByDefault;
    }

    public double getStepSizeThisStep() {
        return stepSizeThisStep;
    }

    public double getStepSizeNextStep() {
        return stepSizeNextStep;
    }

    public double getActivityLevel(int proteinId) {
        return activityLevels.get(proteinId);
    }

    public void resetActivityLevels() {
        Set<Integer> proteinIds = new HashSet<Integer>(activityLevels.keySet());

        for(Integer proteinId : proteinIds) {
            activityLevels.put(proteinId, 0.0);
        }
    }

    //test functions
    private HashMap<Integer, Double> exponentialDecay(HashMap<Integer, Double> quantity) {
        HashMap<Integer, Double> slope = new HashMap<Integer, Double>();
        slope.put(0, (-0.2) * quantity.get(0));
        return slope;
    }

    public void setActivityLevel(int proteinId, double activityLevel) {
        activityLevels.put(proteinId, activityLevel);
    }
}
