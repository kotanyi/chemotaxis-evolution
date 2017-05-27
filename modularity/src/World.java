import org.apache.commons.math3.random.BitsStreamGenerator;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class World implements Cloneable {
    private double input, time, goal, goalClock, averageDesired, rateFactor, trigArg, stochasticGaussianFactor;
    private static double stochasticFoodEquilibrationPeriod;
    private ArrayList<Double> goals, inputs, timesSteadyState, fitnesses, inputsSteadyState, durations, exponentials;
    private ArrayList<ArrayList<Double>> orkunFitnessComponents;
    private ArrayList<ArrayList<Float>> fourierCoefficientsXp, fourierCoefficientsYp;
    private Map<Integer, List<Double>> outputsSteadyState;
    private static List<Double> gammas = new ArrayList<Double>();
    private static List<Gaussian> gaussians;
    private static List<List<Double>> decompositionFunctionsCourses;
    private static double sensitivity, gamma = 1, maxFixationProbability = 0, adaptiveScalingConstant = 15.91, linearScalingConstant = 2.24781, adaptiveDecompositionScalingConstant = 69.7385, adaptiveResponseRangeConstant = 4.38179918544744, orkunLinearWeightingOriginal, orkunAdaptiveWeightingOriginal, fixationProbability;
    private boolean isEquilibrated;
    public static final int DEFAULT_FITNESS = 1;
    private int numberOfConvolutionEstimates;
    private long steps;
    private BitsStreamGenerator random;
    private static final double slope = (Parameters.finalSensitivity - Parameters.initialSensitivity) / (Parameters.numberOfGenerations - Parameters.startIncreasingSensitivityAtGeneration);
    public ArrayList<Double> times, effectorActivities;

    public World(Set<Integer> proteinIds) {
        time = 0;
        goalClock = 0;
        effectorActivities = new ArrayList<Double>();
        times = new ArrayList<Double>();
        outputsSteadyState = new HashMap<Integer, List<Double>>();
        timesSteadyState = new ArrayList<Double>();
        fitnesses = new ArrayList<Double>();
        orkunFitnessComponents = new ArrayList<ArrayList<Double>>();
        decompositionFunctionsCourses = new ArrayList<List<Double>>();
        averageDesired = 0;
        isEquilibrated = false;
        durations = new ArrayList<Double>();
        fourierCoefficientsXp = new ArrayList<ArrayList<Float>>();
        fourierCoefficientsYp = new ArrayList<ArrayList<Float>>();
        rateFactor = 1 - Parameters.stochasticRate * Parameters.stochasticStepSize;
        trigArg = 2 * Math.PI / Parameters.worldLength;
        steps = 0;
        exponentials = new ArrayList<Double>();
        stochasticGaussianFactor = Math.sqrt(2.0 * Parameters.stochasticRate * Parameters.stochasticStepSize * Parameters.stochasticFoodScaling / Parameters.numberOfModes);

        for(int i = 0; i < Parameters.numberOfModes; i++) {
            fourierCoefficientsXp.add(new ArrayList<Float>());
            fourierCoefficientsXp.get(i).add(0f);
            fourierCoefficientsYp.add(new ArrayList<Float>());
            fourierCoefficientsYp.get(i).add(0f);
        }

        for(int i = 0; i < 3; i++) {
            fitnesses.add(0.0);
            orkunFitnessComponents.add(new ArrayList<Double>());
        }

        for(Integer i : proteinIds) {
            outputsSteadyState.put(i, new ArrayList<Double>());
            outputsSteadyState.get(i).add(0.0);
        }

        setInputAndGoal();
        generateDecompositionFunctionsCourses();
    }

    public static void initialiseGaussians() {
        gaussians = new ArrayList<Gaussian>();

        for(int i = 0; i < Parameters.gaussianWidths.length; i++) {
            gaussians.add(new Gaussian(Parameters.gaussianWidths[i], Parameters.gaussianHeights[i], Parameters.gaussianCentres[i]));
        }
    }

    private void setInputAndGoal() {
        if(Parameters.whatToSimulate == 1) {
            if(Parameters.stimulusCourse[2]) {
                input = 0.5;

                if(Parameters.responseFunction[0]) {
                    goal = input;
                } else if(Parameters.responseFunction[1]) {
                    goal = 0.465;
                }
            } else {
                input = 0;

                if(Parameters.responseFunction[0]) {
                    goal = input;
                } else if(Parameters.responseFunction[1]) {
                    goal = 0.5;
                }
            }
        } else {
            if(Parameters.useStochasticFood) {
                input = Parameters.stochasticBaseline;
            } else {
                input = Parameters.gaussianBaselineConcentration;
            }
        }

        inputs = new ArrayList<Double>();
        inputsSteadyState = new ArrayList<Double>();
        inputsSteadyState.add(input);
        goals = new ArrayList<Double>();
    }

    private void generateDecompositionFunctionsCourses() {
        List<Double> estimateTimes = new ArrayList<Double>();

        for(double d = Parameters.stepSize; d <= Parameters.stimulusCourseDuration; d += Parameters.stepSize) {
            estimateTimes.add(d);
        }
        for(int i = 0; i < 3; i++) {
            decompositionFunctionsCourses.add(i, new ArrayList<Double>());

            for(Double d : estimateTimes) {
                time = d;
                advanceInput();

                if(i == 0) {
                    calculateGoalLinear();

                    if(Parameters.stimulusCourse[1]) {
                        goal = goal * linearScalingConstant - linearScalingConstant / 2;
                    } else if(Parameters.stimulusCourse[2]) {
                        goal = ((goal - 0.5) * 2 - 0.137) * 0.8206;
                    }
                }
                else if(i == 1) {
                    calculateGoalAdaptive();

                    if(Parameters.stimulusCourse[1]) {
                        goal = (goal - 0.5) / adaptiveScalingConstant * adaptiveDecompositionScalingConstant;
                    } else if(Parameters.stimulusCourse[2]) {
                        goal = (goal - 0.465) * 13.26 * 0.18407935;
                    }
                } else if(i == 2) {
                    if(Parameters.stimulusCourse[1]) {
                        goal = 1;
                    } else if(Parameters.stimulusCourse[2]) {
                        goal = 1 / Math.sqrt(Parameters.stimulusCourseDuration);
                    }
                }

                decompositionFunctionsCourses.get(i).add(goal);
            }

            time = 0;
            setInputAndGoal();
        }
    }

    @SuppressWarnings("unchecked") //this is to keep java from complaining about the unchecked cast from object to hashmap in w.effectorActivities = (ArrayList<Double>) w.effectorActivities.clone(); (and below)
    public Object clone() throws CloneNotSupportedException {
        World w = null;
        try { w = (World) super.clone(); }
        catch(CloneNotSupportedException e) { System.exit(-47); }
        w.effectorActivities = (ArrayList<Double>) w.effectorActivities.clone();
        w.goals = (ArrayList<Double>) w.goals.clone();
        w.inputs = (ArrayList<Double>) w.inputs.clone();
        w.times = (ArrayList<Double>) w.times.clone();
        w.timesSteadyState = (ArrayList<Double>) w.timesSteadyState.clone();
        w.fitnesses = (ArrayList<Double>) w.fitnesses.clone();
        w.inputsSteadyState = (ArrayList<Double>) w.inputsSteadyState.clone();
        w.orkunFitnessComponents = (ArrayList<ArrayList<Double>>) w.orkunFitnessComponents.clone();
        w.durations = (ArrayList<Double>) w.durations.clone();
        w.exponentials = (ArrayList<Double>) w.exponentials.clone();
        ArrayList<ArrayList<Double>> tmp = new ArrayList<ArrayList<Double>>();

        for(int i = 0; i < 3; i++) {
            ArrayList<Double> j = (ArrayList<Double>) w.orkunFitnessComponents.get(i).clone();
            tmp.add(i, j);
        }
        for(int i = 0; i < 3; i++) {
            w.orkunFitnessComponents.set(i, tmp.get(i));
        }

        w.fourierCoefficientsXp = (ArrayList<ArrayList<Float>>) w.fourierCoefficientsXp.clone();
        w.fourierCoefficientsYp = (ArrayList<ArrayList<Float>>) w.fourierCoefficientsYp.clone();
        ArrayList<ArrayList<Float>> tmpXp = new ArrayList<ArrayList<Float>>();
        ArrayList<ArrayList<Float>> tmpYp = new ArrayList<ArrayList<Float>>();

        for(int i = 0; i < Parameters.numberOfModes; i++) {
            ArrayList<Float> Xp = (ArrayList<Float>) w.fourierCoefficientsXp.get(i).clone();
            tmpXp.add(Xp);
            ArrayList<Float> Yp = (ArrayList<Float>) w.fourierCoefficientsYp.get(i).clone();
            tmpYp.add(Yp);
        }
        for(int i = 0; i < Parameters.numberOfModes; i++) {
            w.fourierCoefficientsXp.set(i, tmpXp.get(i));
            w.fourierCoefficientsYp.set(i, tmpYp.get(i));
        }

        return w;
    }

    public Map<Integer, List<Double>> cloneOutputsSteadyState() {
        Map<Integer, List<Double>> outputsSteadyStateCopy = new HashMap<Integer, List<Double>>();

        for(Integer i : outputsSteadyState.keySet()) {
            outputsSteadyStateCopy.put(i, new ArrayList<Double>());

            for(int j = 0; j < outputsSteadyState.get(i).size(); j++) {
                outputsSteadyStateCopy.get(i).add(outputsSteadyState.get(i).get(j));
            }
        }

        return outputsSteadyStateCopy;
    }

    public void setOutputsSteadyState(Map<Integer, List<Double>> clone) {
        outputsSteadyState = clone;
    }

    //list-related
    public void addEffectorActivity(double effectorActivity) {
        effectorActivities.add(effectorActivity);
    }

    public double getEffectorActivity(int i) {
        return effectorActivities.get(i);
    }

    public double getGoal(int i) {
        return goals.get(i);
    }

    public double getInput(int i) {
        return inputs.get(i);
    }

    public void logTime() {
        times.add(time);
    }

    public double getTime(int i) {
        return times.get(i);
    }

    public int getDataPointCount() {
        return effectorActivities.size();
    }

    public double getOutputsSteadyState(int proteinId, int i) {
        return outputsSteadyState.get(proteinId).get(i);
    }

    public double getInputsSteadyState(int i) {
        return inputsSteadyState.get(i);
    }

    public double getTimesSteadyState(int i) {
        return timesSteadyState.get(i);
    }

    public int getSteadyStateDataPointCount() {
        return timesSteadyState.size();
    }

    public void dropEstimates(Map<Integer, List<Double>> outputs) {
        List<Integer> estimatesToDrop = new ArrayList<Integer>();

        for(int i = 0; i < times.size(); i++) {
            if(times.get(i) % Parameters.stepSize != 0) {
                estimatesToDrop.add(i);
            }
        }

        for(int i = estimatesToDrop.size() - 1; i >= 0; i--) {
            effectorActivities.remove((int) estimatesToDrop.get(i));

            if(Parameters.whatToSimulate == 1 && !Parameters.fitnessCalculation[3]) {
                goals.remove((int) estimatesToDrop.get(i));
            }

            times.remove((int) estimatesToDrop.get(i));
            inputs.remove((int) estimatesToDrop.get(i));

            if(Parameters.whatToSimulate == 1 && Modularity.whatTodo == 3 && !Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
                for(Integer k : outputs.keySet()) {
                    outputs.get(k).remove((int) estimatesToDrop.get(i));
                }
            }
        }
    }

    public static List<Gaussian> getGaussians() {
        return gaussians;
    }

    //sensitivity
    public static void setSensitivity() {
        sensitivity = 1;
    }

    public static double getSensitivity() {
        return sensitivity;
    }

    //input
    public void advanceInput() {
        input = getInput(time);
        inputs.add(input);
    }

    public static double getInput(double time) {
        double input = 0;

        if(Parameters.stimulusCourse[0]) {
            if(time >= 50 && time < 100) {
                input = 1;
            }
            if(time >= 100) {
                input = 0;
            }
        }
        else if(Parameters.stimulusCourse[1]) {
            if(time >= 36 && time <= 86) {
                input = Math.cos((-Math.PI / 2) * ((time - 50 - 36) / (-50))) * Math.cos((-Math.PI / 2) * ((time - 50 - 36) / (-50)));
            }
            else if(time > 86 && time < 156) {
                input = 1;
            }
            else if(time >= 156 && time <= 206) {
                input = Math.cos((Math.PI / 2) * ((time - 50 - 36 - 70) / 50)) * Math.cos((Math.PI / 2) * ((time - 50 - 36 - 70) / 50));
            }
            else if(time > 206) {
                input = 0;
            }
        } else if(Parameters.stimulusCourse[2]) {
            input = Math.sin(time * time) / 2 + 0.5;
        } else if(Parameters.stimulusCourse[3]) {
            if(time > 100 && time <= 175) {
                input = Math.cos(-Math.PI / 2 * (time - 175) / -74) * Math.cos(-Math.PI / 2 * (time - 175) / -74) / 2;
            } else if(time > 175 && time <= 250) {
                input = Math.cos(Math.PI / 2 * (time - 175 - 1) / 74) * Math.cos(Math.PI / 2 * (time - 175 - 1) / 74) / 2;
            } else if(time > 250 && time <= 350) {
                input = 0;
            } else if(time > 350 && time <= 400) {
                input = Math.cos(-Math.PI / 4 * (time - 175 - 176 - 49) / -49) * Math.cos(-Math.PI / 4 * (time - 175 - 176 - 49) / -49) - 0.5;

                if(input < 0) {
                    input = 0;
                }
            } else if(time > 400 && time <= 625) {
                input = 0.5;
            } else if(time > 625 && time <= 700) {
                input = Math.cos(Math.PI / 2 * (time - 175 - 451) / 74) * Math.cos(Math.PI / 2 * (time - 175 - 451) / 74) / 3 + ((double) 1/6);
            } else if(time > 700 && time <= 775) {
                input = Math.cos(-Math.PI / 2 * (time - 175 - 526 - 74) / -74) * Math.cos(-Math.PI / 2 * (time - 175 - 526 - 74) / -74) / 3 + ((double) 1/6);
            } else if(time > 775) {
                input = 0.5;
            }
        } else if(Parameters.stimulusCourse[4]) {
            if(time > 1100 && time <= 1175) {
                input = Math.cos(-Math.PI / 2 * (time - 1175) / -75) * Math.cos(-Math.PI / 2 * (time - 1175) / -75) / 2;
            } else if(time > 1175 && time <= 1250) {
                input = Math.cos(Math.PI / 2 * (time - 1175) / 75) * Math.cos(Math.PI / 2 * (time - 1175) / 75) / 2;
            } else if(time > 1250 && time <= 1350) {
                input = 0;
            } else if(time > 1350 && time <= 1400) {
                input = Math.cos(-Math.PI / 4 * (time - 1400) / -50) * Math.cos(-Math.PI / 4 * (time - 1400) / -50) - 0.5;
            } else if(time > 1400 && time <= 1750) {
                input = 0.5;
            } else if(time > 1750 && time <= 1825) {
                input = Math.cos(Math.PI / 2 * (time - 1750) / 75) * Math.cos(Math.PI / 2 * (time - 1750) / 75) / 3 + ((double) 1/6);
            } else if(time > 1825 && time <= 1900) {
                input = Math.cos(-Math.PI / 2 * (time - 1900) / -75) * Math.cos(-Math.PI / 2 * (time - 1900) / -75) / 3 + ((double) 1/6);
            } else if(time > 1900) {
                input = 0.5;
            }
        } else if(Parameters.stimulusCourse[5]) {
            if(time == 1) {
                input = 1;
            } else {
                input = 0;
            }
        } else if(Parameters.stimulusCourse[6]) {
            input = Math.exp(-(time - Parameters.stimulusCourseDuration / 2) * (time - Parameters.stimulusCourseDuration / 2) / (2 * (Parameters.stimulusCourseDuration / 10) * (Parameters.stimulusCourseDuration / 10)));
        } else if(Parameters.stimulusCourse[7]) {
            if(time >= 191.71 && time <= 1176.1) {
                input = 1;
            } else {
                input = 0;
            }
        }

        return input;
    }

    public double getInput() {
        return input;
    }

    //goal
    public void calculateGoal() {
        if(Parameters.responseFunction[0]) {
            calculateGoalLinear();
        }
        else if(Parameters.responseFunction[1]) {
            calculateGoalAdaptive();
        }
    }

    private void calculateGoalLinear() {
        goal = input;
        goals.add(goal);
    }

    private void calculateGoalAdaptive() {
        if(Parameters.stimulusCourse[0]) {
            if(time >= 50 && time < 100)
                goal = 0.5 - (Math.exp(-(time - 50) / Parameters.timeConstant) / 2);
            if(time >= 100)
                goal = 0.5 + (Math.exp(-(time - 100) / Parameters.timeConstant) / 2);
        }
        else if(Parameters.stimulusCourse[1]) {
            if(time < 36) {
                goal = 0.5;
            }
            else if(time >= 36 && time <= 86) {
                goal = -(-Math.PI / 50) * Math.cos((-Math.PI / 2) * ((time - 50 - 36) / (-50))) * Math.sin((-Math.PI / 2) * ((time - 50 - 36) / (-50))) * adaptiveScalingConstant + 0.5;
            }
            else if(time > 86 && time < 156) {
                goal = 0.5;
            }
            else if(time >= 156 && time <= 206) {
                goal = -(-Math.PI / 50) * Math.cos((Math.PI / 2) * ((time - 50 - 36 - 70) / 50)) * Math.sin((Math.PI / 2) * ((time - 50 - 36 - 70) / 50)) * adaptiveScalingConstant + 0.5;
            }
            else if(time > 206) {
                goal = 0.5;
            }
        } else if(Parameters.stimulusCourse[2]) {
            goal = 2 * time * Math.cos(time * time) / 13.26 + 0.465;
        }

        goals.add(goal);
    }

    //fitness calculation
    public void calculateFitness() {
        if(Parameters.whatToSimulate == 1) {
            if(Parameters.fitnessCalculation[0]) {
                leastSquares();
            } else if(Parameters.fitnessCalculation[1]) {
                rewardingForFit();
            } else if(Parameters.fitnessCalculation[2]) {
                decompose();
            } else if(Parameters.fitnessCalculation[3]) {
                orkun();
            }
        } else if(Parameters.whatToSimulate == 2) {
            if(Parameters.curveFitting) {
                leastSquares2();
            } else {
                averageHarvestDuration();
            }
        }
    }

    private void leastSquares() {
        double fitness = 0;

        for(int i = 0; i < goals.size(); i++)
            fitness += (effectorActivities.get(i) - goals.get(i)) * (effectorActivities.get(i) - goals.get(i));

        //fitness = 1/fitness;
        fitness = 1 - fitness/averageDesired;
        fitnesses.set(DEFAULT_FITNESS, fitness);
    }

    private void leastSquares2() {
        int startFrom = 0;
        double sum = 0, activity;

        for(int i = 0; i < Modularity.packer.get(0).size(); i++) {
            for(int j = startFrom; j < times.size()-1; j++) {
                if(Modularity.packer.get(0).get(i) >= times.get(j) && Modularity.packer.get(0).get(i) < times.get(j+1)) {
                    activity = effectorActivities.get(j) + (effectorActivities.get(j+1) - effectorActivities.get(j)) / (times.get(j+1) - times.get(j)) * (Modularity.packer.get(0).get(i) - times.get(j));
                    sum += (activity - Modularity.packer.get(1).get(i)) * (activity - Modularity.packer.get(1).get(i));
                    startFrom = j;
                    break;
                }
            }
        }

        fitnesses.set(DEFAULT_FITNESS, 1 / sum);
    }

    private void rewardingForFit() {
        double fitness = 0;

        for(int i = 0; i < goals.size(); i++) {
            fitness += Math.exp(-(effectorActivities.get(i) - goals.get(i)) * (effectorActivities.get(i) - goals.get(i)) / (2 * Parameters.gaussianFunctionParameterC * Parameters.gaussianFunctionParameterC));
        }

        fitness = 1 - averageDesired/fitness;
        fitnesses.set(DEFAULT_FITNESS, fitness);
    }

    private void decompose() {
        List<List<Double>> effectorActivitiesScaled = new ArrayList<List<Double>>();

        for(int i = 0; i < 3; i++) {
            effectorActivitiesScaled.add(i, new ArrayList<Double>());
        }

        for(Double d : effectorActivities) {
            if(Parameters.stimulusCourse[1]) {
                effectorActivitiesScaled.get(0).add(d * linearScalingConstant);
                effectorActivitiesScaled.get(1).add(d * adaptiveResponseRangeConstant);
            } else if(Parameters.stimulusCourse[2]) {
                effectorActivitiesScaled.get(0).add(d * 2 * 0.8206);
                effectorActivitiesScaled.get(1).add(d * 13.26 * 0.18407935);
            }

            effectorActivitiesScaled.get(2).add(d);
        }

        for(int i = 0; i < effectorActivitiesScaled.size(); i++) {
            for(int j = 0; j < effectorActivitiesScaled.get(i).size(); j++) {
                if(Parameters.stimulusCourse[1]) {
                    fitnesses.set(i, fitnesses.get(i) + effectorActivitiesScaled.get(i).get(j) * decompositionFunctionsCourses.get(i).get(j) / (Parameters.stimulusCourseDuration / Parameters.stepSize));
                } else if(Parameters.stimulusCourse[2]) {
                    fitnesses.set(i, fitnesses.get(i) + effectorActivitiesScaled.get(i).get(j) * decompositionFunctionsCourses.get(i).get(j) * 0.015);
                }
            }
        }
    }

    private void orkun() {
        double deltaF, averageTumbling, optimalTumbling = 0.45, fitness = 0, averageF;

        for(int i = 1; i < getDataPointCount(); i++) {
            deltaF = inputs.get(i) - inputs.get(i-1);
            averageTumbling = (effectorActivities.get(i) + effectorActivities.get(i-1)) / 2;
            averageF = (inputs.get(i) + inputs.get(i-1)) / 2;
            orkunFitnessComponents.get(0).add(Parameters.orkunLinearWeighting * averageF * (averageTumbling - optimalTumbling));
            orkunFitnessComponents.get(1).add(Parameters.orkunAdaptiveWeighting * deltaF * (optimalTumbling - averageTumbling));
            orkunFitnessComponents.get(2).add((optimalTumbling - averageTumbling) * (optimalTumbling - averageTumbling));
            fitness += Parameters.orkunAdaptiveWeighting * deltaF * (optimalTumbling - averageTumbling) + Parameters.orkunLinearWeighting * averageF * (averageTumbling - optimalTumbling) - (optimalTumbling - averageTumbling) * (optimalTumbling - averageTumbling);
        }

        fitnesses.set(DEFAULT_FITNESS, fitness);
    }

    public static void transitionFromLinearToAdaptive(int generation) {
        if(generation == 1) {
            orkunLinearWeightingOriginal = Parameters.orkunLinearWeighting;
            orkunAdaptiveWeightingOriginal = Parameters.orkunAdaptiveWeighting;
            Parameters.orkunAdaptiveWeighting = 0;
        } else if(generation >= Parameters.orkunTransitionAt && generation < Parameters.orkunTransitionAt + Parameters.orkunTransitionDuration) {
            Parameters.orkunLinearWeighting = orkunLinearWeightingOriginal * (1 - (generation - Parameters.orkunTransitionAt) / (double) Parameters.orkunTransitionDuration);
            Parameters.orkunAdaptiveWeighting = orkunAdaptiveWeightingOriginal * (generation - Parameters.orkunTransitionAt) / (double) Parameters.orkunTransitionDuration;
        } else if(generation == Parameters.orkunTransitionAt + Parameters.orkunTransitionDuration) {
            Parameters.orkunLinearWeighting = 0;
            Parameters.orkunAdaptiveWeighting = orkunAdaptiveWeightingOriginal;
        }
    }

    public void setFitness(double fitness) {
        fitnesses.set(DEFAULT_FITNESS, fitness);
    }

    public double getFitness() {
        return fitnesses.get(World.DEFAULT_FITNESS);
    }

    public double getFitness(int i) {
        return fitnesses.get(i);
    }

    public void calculateAverageDesired() {
        if(Parameters.fitnessCalculation[0]) {
            for(int i = 0; i < goals.size(); i++)
                averageDesired += (goals.get(i) - Parameters.baselineFitness) * (goals.get(i) - Parameters.baselineFitness);
        }
        else if(Parameters.fitnessCalculation[1]) {
            for(int i = 0; i < goals.size(); i++)
                averageDesired += Math.exp(-(goals.get(i) - Parameters.baselineFitness) * (goals.get(i) - Parameters.baselineFitness) / (2 * Parameters.gaussianFunctionParameterC * Parameters.gaussianFunctionParameterC));
        }
    }

    public double getOrkunFitnessComponents(int i, int j) {
        return orkunFitnessComponents.get(i).get(j);
    }

    //fitness functions
    public static int fitnessFunction(double goodnessOfFitPre, double goodnessOfFitPost) {
        int acceptanceProbability = -1;

        if(Parameters.fitnessFunction[0]) {
            acceptanceProbability = hillClimbing(goodnessOfFitPre, goodnessOfFitPost);
        } else if(Parameters.fitnessFunction[1] || Parameters.fitnessFunction[2]) {
            acceptanceProbability = kimura(goodnessOfFitPre, goodnessOfFitPost);
        //previously, non-strict was automatically used for decomposition and orkun and strict was used for least squares and rewarding for fit. now it depends on the user-set parameters file.
        } else if(Parameters.fitnessFunction[3]) {
            acceptanceProbability = simulatedAnnealing(goodnessOfFitPre, goodnessOfFitPost);
        } else if(Parameters.fitnessFunction[4]) {
            acceptanceProbability = simulatedAnnealingStrict(goodnessOfFitPre, goodnessOfFitPost);
        }

        return acceptanceProbability;
    }

    private static int kimura(double goodnessOfFitPre, double goodnessOfFitPost) {
        double fixationProbability, random = Modularity.random.nextDouble();

        if(goodnessOfFitPost - goodnessOfFitPre == 0 || goodnessOfFitPre == 0) {
            fixationProbability = 1 / ((double) Parameters.populationSize);
        }
        else {
            double s = (goodnessOfFitPost - goodnessOfFitPre) / goodnessOfFitPre;
            fixationProbability = (1 - Math.exp(-2 * s)) / (1 - Math.exp(-2 * ((double) Parameters.populationSize) * s));
        }

        if(Parameters.fitnessFunction[2]) {
            fixationProbability = scale(fixationProbability);
        }

        if(random <= fixationProbability) {
            return 1;
        }
        else {
            return 0;
        }
    }

    private static double scale(double fixationProbability) {
        if(fixationProbability > maxFixationProbability) {
            maxFixationProbability = fixationProbability;
        }

        fixationProbability *= gamma;

        if(fixationProbability > 1) {
            System.out.println("fixationProbability is greater than 1!");
        }

        return fixationProbability;
    }

    public static void setGamma(int generation) {
        if(generation > 0 && generation % 1000 == 0) {
            gammas.add(gamma);
            gamma = -(Math.ceil(Math.log10(maxFixationProbability)) + 1);

            if(gamma == -1) {
                gamma = 1;
            }
            else {
                gamma = Math.pow(10, gamma);
            }

            maxFixationProbability = 0;
        }
    }

    private static int simulatedAnnealingStrict(double goodnessOfFitPre, double goodnessOfFitPost) {
        if(goodnessOfFitPost - goodnessOfFitPre >= 0) {
            fixationProbability = 1;
            return 1;
        } else if(goodnessOfFitPost - goodnessOfFitPre < 0 && goodnessOfFitPost >= 0) {
            double random = Modularity.random.nextDouble();

            if(goodnessOfFitPre == 0) { //if the denominator in the selection coefficient expression is equal to zero, replace with a small number
                goodnessOfFitPre = 1e-6;
            }

            fixationProbability = Math.exp(((goodnessOfFitPost - goodnessOfFitPre) / goodnessOfFitPre) / Parameters.simulatedAnnealingTemperature);

            if(random <= fixationProbability) {
                return 1;
            } else {
                return 0;
            }
        } else {
            fixationProbability = 0;
            return 0;
        }
    }

    private static int simulatedAnnealing(double goodnessOfFitPre, double goodnessOfFitPost) {
        if(goodnessOfFitPost - goodnessOfFitPre >= 0) {
            fixationProbability = 1;
            return 1;
        } else {
            double random = Modularity.random.nextDouble();

            if(goodnessOfFitPre == 0) {
                goodnessOfFitPre = 1e-6;
            }

            fixationProbability = Math.exp(((goodnessOfFitPost - goodnessOfFitPre) / goodnessOfFitPre) / Parameters.simulatedAnnealingTemperature);

            if(random <= fixationProbability) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public static int simulatedAnnealingFlat(double goodnessOfFitPre, double goodnessOfFitPost, List<Double> otherFitnesses) { //otherFitnesses[0] = linear pre, [1] = linear post, [2] = adaptive pre, [3] = adaptive post
        boolean otherResponsesAreWeak = Math.abs(otherFitnesses.get(0)) < 0.07 && Math.abs(otherFitnesses.get(1)) < 0.07;

        if(Math.abs(goodnessOfFitPost - 0.5) <= Math.abs(goodnessOfFitPre - 0.5) && otherResponsesAreWeak) {
            return 1;
        } else if(Math.abs(goodnessOfFitPost - 0.5) > Math.abs(goodnessOfFitPre - 0.5) && otherResponsesAreWeak) {
            double random = Modularity.random.nextDouble();

            if(goodnessOfFitPre == 0.5) { //this is necessary to ensure that Math.abs(goodnessOfFitPre - 0.5) is not zero - see below
                goodnessOfFitPre += 1e-6;
            }

            if(random <= Math.exp(((Math.abs(goodnessOfFitPre - 0.5) - Math.abs(goodnessOfFitPost - 0.5)) / Math.abs(goodnessOfFitPre - 0.5)) / Parameters.simulatedAnnealingTemperature)) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    public static int simulatedAnnealingFlatAdaptive(double goodnessOfFitPre, double goodnessOfFitPost, List<Double> otherFitnesses) {
        boolean otherResponsesAreWeak = Math.abs(otherFitnesses.get(0)) < 0.2 && Math.abs(otherFitnesses.get(1)) < 0.2 && goodnessOfFitPre > 0.48 && goodnessOfFitPre < 0.52 && goodnessOfFitPost > 0.48 && goodnessOfFitPost < 0.52;

        if(otherFitnesses.get(3) - otherFitnesses.get(2) >= 0 && otherResponsesAreWeak) {
            return 1;
        } else if(otherFitnesses.get(3) - otherFitnesses.get(2) < 0 && otherFitnesses.get(3) >= 0 && otherResponsesAreWeak) {
            double random = Modularity.random.nextDouble();

            if(otherFitnesses.get(2) == 0) {
                otherFitnesses.set(2, 1e-6);
            }

            if(random <= Math.exp(((otherFitnesses.get(3) - otherFitnesses.get(2)) / otherFitnesses.get(2)) / Parameters.simulatedAnnealingTemperature)) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    private static int hillClimbing(double goodnessOfFitPre, double goodnessOfFitPost) {
        if(goodnessOfFitPost - goodnessOfFitPre >= 0) {
            fixationProbability = 1;
            return 1;
        }
        else {
            fixationProbability = 0;
            return 0;
        }
    }

    //steady state
    public void waitUntilSteadyState(RungeKutta rk, ProteinNetwork pn, Bug bug) {
        double toBeAddedToActivitiesOlder;
        boolean steadyState = true;
        Map<Integer, Double> stdDeviationsRecent = new HashMap<Integer, Double>(), stdDeviationsOlder = new HashMap<Integer, Double>();
        List<Double> times = new ArrayList<Double>();
        times.add(time);
        Map<Integer, Queue<Double>> activitiesRecent = new HashMap<Integer, Queue<Double>>(), activitiesOlder = new HashMap<Integer, Queue<Double>>();

        for(Integer i : pn.getProteinIds()) {
            activitiesRecent.put(i, new LinkedList<Double>());
            activitiesRecent.get(i).offer(0.0);
            activitiesOlder.put(i, new LinkedList<Double>());
        }

        while(true) {
            if(time >= Parameters.steadyStateMaxEstimates) {
                if(Modularity.whatTodo == 3) {
                    System.out.println("Steady state has not been reached, proceeding to the stimulus course.");
                }

                break;
            }
            if(!rk.cashKarp(pn, this, bug, true)) {
                System.out.println("waitUntilSteadyState: fourthOrder used!");
                System.exit(-47);
            }

            time += rk.getStepSizeThisStep();

            if(time % Parameters.stepSize == 0) {
                times.add(time);
            } else {
                continue;
            }

            inputsSteadyState.add(input);

            for(Integer i : pn.getProteinIds()) {
                outputsSteadyState.get(i).add(rk.getActivityLevel(i));

                if(activitiesRecent.get(i).size() == Parameters.steadyStateQueueSize) {
                    toBeAddedToActivitiesOlder = activitiesRecent.get(i).poll();
                    activitiesRecent.get(i).offer(rk.getActivityLevel(i));

                    if(time % 20 == 0) {
                        stdDeviationsRecent.put(i, fastStdDeviation(activitiesRecent.get(i)));
                    }

                    //the if clause with < has to go after the one with == otherwise rk.getActivityLevel(i) would be added (when activitiesOlder.get(i).size() == 9) and subsequently the clause with == would be executed resulting in subtracting and adding rk.getActivityLevel(i) again
                    if(activitiesOlder.get(i).size() == Parameters.steadyStateQueueSize) {
                        activitiesOlder.get(i).poll();
                        activitiesOlder.get(i).offer(toBeAddedToActivitiesOlder);

                        if(time % 20 == 0) {
                            stdDeviationsOlder.put(i, fastStdDeviation(activitiesOlder.get(i)));
                        }
                    }
                    if(activitiesOlder.get(i).size() < Parameters.steadyStateQueueSize) {
                        activitiesOlder.get(i).offer(toBeAddedToActivitiesOlder);

                        if(activitiesOlder.get(i).size() == Parameters.steadyStateQueueSize) {
                            stdDeviationsOlder.put(i, fastStdDeviation(activitiesOlder.get(i)));
                        }
                    }
                }
                if(activitiesRecent.get(i).size() < Parameters.steadyStateQueueSize) {
                    activitiesRecent.get(i).offer(rk.getActivityLevel(i));

                    if(activitiesRecent.get(i).size() == Parameters.steadyStateQueueSize) {
                        stdDeviationsRecent.put(i, fastStdDeviation(activitiesRecent.get(i)));
                    }
                }
            }

            if(!stdDeviationsOlder.isEmpty()) {
                for(Integer i : pn.getProteinIds()) {
                    if(!(Math.abs(stdDeviationsOlder.get(i) - stdDeviationsRecent.get(i)) <= Parameters.steadyStateThreshold)) {
                        steadyState = false;
                    }
                }
                if(steadyState) {
                    isEquilibrated = true;
                    break;
                } else {
                    steadyState = true;
                }
            }
        }

        for(int i = 0; i < times.size(); i++) {
            timesSteadyState.add(i, -(times.get(times.size()-1) - times.get(i)));
        }

        time = 0;
        rk.resetStepSizeNextStep();
    }

    private double calculateStdDeviation(Queue<Double> activities) {
        double sum = 0, average, squaredSum = 0, stdDeviation;

        for(Double d : activities) {
            sum += d;
        }

        average = sum/activities.size();

        for(Double d : activities) {
            squaredSum += (d - average) * (d - average);
        }

        stdDeviation = Math.sqrt(squaredSum / activities.size());
        return stdDeviation;
    }

    public double fastStdDeviation(Collection<Double> activities) {
        int k = 1;
        double Mk = 0, MkPrevious, Qk = 0, stdDeviation;

        for(Double x : activities) {
            if(k == 1) {
                Mk = x;
                Qk = 0;
            } else {
                MkPrevious = Mk;
                Mk = Mk + (x - Mk) / k;
                Qk = Qk + ((k - 1) * (x - MkPrevious) * (x - MkPrevious)) / k;
            }

            k++;
        }

        stdDeviation = Math.sqrt(Qk / activities.size());
        return stdDeviation;
    }

    public boolean isEquilibrated() {
        return isEquilibrated;
    }

    //time
    public void increaseTime(double increment) {
        time += increment;
        goalClock += increment;

        if(time - (int) time == 0 && goalClock - (int) goalClock != 0) {
            goalClock = Math.round(goalClock);
        }
    }

    public void setTime(double time) {
        this.time = time;
    }

    public double getTime() {
        return time;
    }

    //print
    public void appendToFile(int generation, boolean networkShouldBePrinted, double fitnessRejected) {
        NumberFormat format = new DecimalFormat("0.######E0");

        if(networkShouldBePrinted) {
            if(Parameters.fitnessCalculation[2]) {
                Modularity.fitnesses.println(generation + "," + format.format(this.fitnesses.get(0)) + "," + format.format(this.fitnesses.get(1)) + "," + format.format(this.fitnesses.get(2)));
            } else {
                Modularity.fitnesses.println(generation + "," + format.format(this.fitnesses.get(DEFAULT_FITNESS)));
            }

            Modularity.fitnessRejectedWriter.println(generation + "," + format.format(fitnessRejected));
            Modularity.fixationProbabilities.println(generation + "," + format.format(fixationProbability));

            if(!Parameters.printPreAndPostDurations && Parameters.whatToSimulate == 2) {
                appendToDurations(generation);
            }
        } if(outputShouldBePrinted(generation) && Parameters.whatToSimulate != 2) {
            for(int i = 0; i < effectorActivities.size(); i++) {
                Modularity.output.println(generation + "," + format.format(times.get(i)) + "," + format.format(effectorActivities.get(i)));
            }

            if(Parameters.fitnessCalculation[3]) {
                if(isEquilibrated) {
                    for(int i = 0; i < orkunFitnessComponents.get(0).size(); i++) {
                        Modularity.orkunComponents.println(generation + "," + format.format(times.get(i+1)) + "," + format.format(orkunFitnessComponents.get(0).get(i)) + "," + format.format(orkunFitnessComponents.get(1).get(i)) + "," + format.format(orkunFitnessComponents.get(2).get(i)));
                    }
                } else {
                    Modularity.orkunComponents.println(generation + ",\"The network that was accepted in this generation has not reached an equilibrium and the stimulus course was therefore not applied.\"");
                }
            }
        }
    }

    public void appendToDurations(int generation) {
        Modularity.durations.print(generation);

        for(int i = 0; i < durations.size(); i++) {
            Modularity.durations.print("," + durations.get(i).intValue());
        }

        Modularity.durations.println();
    }

    private boolean outputShouldBePrinted(int generation) {
        if(generation == Parameters.numberOfGenerations) {
            return true;
        } else if(Parameters.orkunTransition && generation == Parameters.orkunTransitionAt - 1) {
            return true;
        } else {
            return false;
        }
    }

    //bug
    public boolean isTimeWhole(int howManyStepsInUnitTime) {
        if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            return steps % howManyStepsInUnitTime == 0;
        } else {
            return time % Parameters.stepSize == 0;
        }
    }

    public int getPreviousWholeTimeIndex() {
        int i = times.size() - 2, j;

        for(j = i; j >= 0; j--) {
            if(times.get(j) % Parameters.stepSize == 0) {
                break;
            }
        }

        if(j == -1)
            j = 0;

        return j;
    }

    public List<List<Double>> getPartOfList(int indexStart) {
        List<List<Double>> listPart = new ArrayList<List<Double>>();
        listPart.add(new ArrayList<Double>());
        listPart.add(new ArrayList<Double>());

        for(int i = indexStart; i < times.size(); i++) {
            listPart.get(0).add(times.get(i));
            listPart.get(1).add(effectorActivities.get(i));
        }

        return listPart;
    }

    public static List<Double> getInputGaussians(double position, boolean isTumbling, double velocity, double driftVelocity, double stepSize) {
        List<Double> sums = new ArrayList<Double>(), positions = new ArrayList<Double>();

        if(isTumbling) {
            velocity = driftVelocity;
        } else {
            velocity += driftVelocity;
        }

        positions.add(position);
        positions.add(position + velocity * stepSize / 5);
        positions.add(position + velocity * stepSize * 3 / 10);
        positions.add(position + velocity * stepSize * 3 / 5);
        positions.add(position + velocity * stepSize);
        positions.add(position + velocity * stepSize * 7 / 8);

        for(int i = 0; i < 6; i++) {
            sums.add(i, Parameters.gaussianBaselineConcentration);

            for(Gaussian g : gaussians) {
                sums.set(i, sums.get(i) + g.getFoodAt(positions.get(i)));
            }
        }

        return sums;
    }

    public double getFoodAt(double position) {
        double food = Parameters.gaussianBaselineConcentration;

        for(Gaussian g : gaussians) {
            food += g.getFoodAt(position);
        }

        return food;
    }

    public void logFood(double position) {
        inputs.add(getFoodAt(position));
    }

    private void averageHarvestDuration() {
        double averageDuration = 0;

        for(double d : durations) {
            averageDuration += 1/d;
        }

        averageDuration /= Parameters.numberOfConcurrentBugs;
        fitnesses.set(DEFAULT_FITNESS, averageDuration);
    }

    public void setDurations(ArrayList<Double> durations) {
        this.durations = durations;
    }

    public double stochasticFoodDistribution(double position, int index) {
        double sumXp = 0, sumYp = 0, trigArgPosition = trigArg * position, trigArgComplete, sum;

        if(index == -1) {
            index = fourierCoefficientsXp.get(0).size() - 1;
        } else {
            addFourierCoefficientsIfNecessary(index);
        }

        for(int i = 0; i < Parameters.numberOfModes; i++) {
            trigArgComplete = trigArgPosition * (i + 1);
            sumXp += fourierCoefficientsXp.get(i).get(index) * Math.cos(trigArgComplete);
            sumYp += fourierCoefficientsYp.get(i).get(index) * Math.sin(trigArgComplete);
        }

        sum = Parameters.stochasticBaseline + sumXp + sumYp;

        if(Parameters.stochasticIsFoodSquared) {
            return sum * sum;
        } else {
            return sum > 0 ? sum : 0;
        }
    }

    public void nextFourierCoefficient() {
        int lastEntryId = fourierCoefficientsXp.get(0).size() - 1;

        for(int i = 0; i < Parameters.numberOfModes; i++) {
            fourierCoefficientsXp.get(i).add((float) (fourierCoefficientsXp.get(i).get(lastEntryId) * rateFactor + random.nextGaussian() * stochasticGaussianFactor));
            fourierCoefficientsYp.get(i).add((float) (fourierCoefficientsYp.get(i).get(lastEntryId) * rateFactor + random.nextGaussian() * stochasticGaussianFactor));
        }
    }

    public double logStochasticFood(double position) {
        double food = stochasticFoodDistribution(position, World.getIndexOfNearestStochasticTime(time));
        inputs.add(food);
        return food;
    }

    public boolean addFourierCoefficientsIfNecessary(int id) {
        int lastEntryId = fourierCoefficientsXp.get(0).size() - 1;

        for(int i = 0; i < id - lastEntryId; i++) {
            nextFourierCoefficient();
        }

        return (id - lastEntryId > 0);
    }

    public List<ArrayList<ArrayList<Float>>> getFourierCoefficients() {
        List<ArrayList<ArrayList<Float>>> fourierCoefficients = new ArrayList<ArrayList<ArrayList<Float>>>();
        fourierCoefficients.add(fourierCoefficientsXp);
        fourierCoefficients.add(fourierCoefficientsYp);
        return fourierCoefficients;
    }

    public void setFourierCoefficients(List<ArrayList<ArrayList<Float>>> fourierCoefficients) {
        fourierCoefficientsXp = fourierCoefficients.get(0);
        fourierCoefficientsYp = fourierCoefficients.get(1);
    }

    public static int getIndexOfNearestStochasticTime(double time) {
        if(Parameters.equilibrateStochasticFood)
            time += World.stochasticFoodEquilibrationPeriod;

        double d = time / Parameters.stochasticStepSize;
        int floor = (int) d;

        if(d < floor + 0.5) {
            return floor;
        } else {
            return floor + 1;
        }
    }

    public void equilibrateStochasticFood() {
        if(!Parameters.averageStochasticFood) {
            for(double d = 0; d <= World.stochasticFoodEquilibrationPeriod; d += Parameters.stochasticStepSize) {
                nextFourierCoefficient();
            }
        } else {
            for(double d = 0; d <= World.stochasticFoodEquilibrationPeriod; d += Parameters.stochasticStepSize) {
                fourierCoefficientsXp.get(0).add(0f);
                fourierCoefficientsYp.get(0).add((float) (6.5 / 100 * 2 * Math.PI));
            }
        }
    }

    public static void determineEquilibrationPeriod(double realisticBugMemoryLengthPre, double realisticBugMemoryLengthPost) {
        List<Double> l = new ArrayList<Double>();
        l.add(1 / Parameters.stochasticRate);
        l.add(realisticBugMemoryLengthPre * Parameters.memoryLengthFactor);
        l.add(realisticBugMemoryLengthPost * Parameters.memoryLengthFactor);
        World.stochasticFoodEquilibrationPeriod = Collections.max(l);
    }

    public Queue<Double> getFoodsPast(ProteinNetwork pn, Bug bug) {
        Queue<Double> foodsPast = new LinkedList<Double>();

        for(int i = -numberOfConvolutionEstimates + 1; i <= 0; i++) {
            if((!Parameters.equilibrateStochasticFood || Parameters.whatToSimulate == 1) || Parameters.curveFitting) {
                foodsPast.offer(0.0);
            } else {
                foodsPast.offer(stochasticFoodDistribution(bug.getPosition(), getIndexOfNearestStochasticTime(i * pn.bugStepSize)));
            }
        }

        return foodsPast;
    }

    public void increaseSteps() {
        if(Parameters.useRealisticBug) {
            steps += 2;
        } else {
            steps++;
        }
    }

    public long getSteps() {
        return steps;
    }

    public void calculateExponentials(ProteinNetwork pn) {
        double time, estimate, tau = pn.realisticBugParameters.get(3);
        double aScaled = pn.realisticBugParameters.get(0) / tau;
        double bScaled = pn.realisticBugParameters.get(1) / (tau * tau);
        double cScaled = pn.realisticBugParameters.get(2) / (tau * tau * tau);

        for(int i = numberOfConvolutionEstimates - 1; i >= 0; i--) {
            time = i * pn.bugStepSize;
            estimate = (aScaled + bScaled * time + cScaled * time * time) * Math.exp(-time / tau);

            if(i != numberOfConvolutionEstimates - 1 && i != 0 ) {
                if(i % 2 == 1) {
                    estimate *= 4;
                } else {
                    estimate *= 2;
                }
            }

            exponentials.add(estimate);
        }
    }

    public List<Double> getExponentials() {
        return exponentials;
    }

    public void resetExponentials() {
        exponentials.clear();
    }

    public void setRandom(BitsStreamGenerator r) {
        random = r;
    }

    public void addInput(double input) {
        inputs.add(input);
    }

    public void memoryCalculations(ProteinNetwork pn) {
        numberOfConvolutionEstimates = (int) Math.floor(pn.realisticBugParameters.get(3) * Parameters.memoryLengthFactor / pn.bugStepSize);

        if(numberOfConvolutionEstimates % 2 == 0) {
            numberOfConvolutionEstimates--;
        }
    }

    public static void updateSensitivity(int generation) {
        if(Parameters.constrainSensitivity && generation >= Parameters.startIncreasingSensitivityAtGeneration) {
            Parameters.sensitivity = slope * (generation - Parameters.startIncreasingSensitivityAtGeneration) + Parameters.initialSensitivity;
        } else if(Parameters.constrainSensitivityInvertedEncodedByA && generation >= Parameters.startIncreasingSensitivityAtGeneration) {
            Parameters.mutationRateOriginal[11] = 5;
            Parameters.minimumConstantAllowed = -Double.MAX_VALUE;
            Parameters.minimumLinearAllowed = -Double.MAX_VALUE;
            Parameters.sensitivity = slope * (generation - Parameters.startIncreasingSensitivityAtGeneration) + Parameters.initialSensitivity;
        }
    }

    public void printFourierCoefficients() {
        PrintWriter p;

        try {
            p = new PrintWriter(new BufferedWriter(new FileWriter("fourier_coefficients")));

            for(int i = 0; i < fourierCoefficientsXp.get(0).size(); i++) {
                for(int j = 0; j < Parameters.numberOfModes; j++) {
                    p.print(fourierCoefficientsXp.get(j).get(i) + ",");
                }
                for(int j = 0; j < Parameters.numberOfModes; j++) {
                    if(j != Parameters.numberOfModes-1) {
                        p.print(fourierCoefficientsYp.get(j).get(i) + ",");
                    } else {
                        p.println(fourierCoefficientsYp.get(j).get(i));
                    }
                }
            }

            p.close();
        } catch(IOException e) {}
    }
}
