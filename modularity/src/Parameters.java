import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

public class Parameters {
    public static int initialProteinCount, maxNumberOfIterations, populationSize, maxNumberOfProteins, numberOfGenerations, deletionToAdditionRatio, minNumberOfProteins, orkunTransitionAt, orkunTransitionDuration, steadyStateQueueSize, steadyStateMaxEstimates, printEvery, whatToSimulate, numberOfConcurrentBugs, howManyBugsToConsider, maxWait, foodRunningSumLimit, printGenerationEvery, loadNetwork, flushEvery, numberOfModes, startIncreasingSensitivityAtGeneration;
    public static double[] mutationRateOriginal, mutationRate, gaussianWidths, gaussianHeights, gaussianCentres;
    public static double advanceTimeByDefault, absoluteTolerance, relativeTolerance, timeConstant, advanceTimeByMax, gaussianFunctionParameterC, maxInteractionStrength, baselineFitness, stimulusCourseDuration, simulatedAnnealingTemperature, orkunAdaptiveWeighting, orkunLinearWeighting, steadyStateThreshold, worldLength, effectorMotorAssociationRate, effectorMotorDissociationRate, stepSize, gaussianDetectionThreshold, maxDurationFactor, driftVelocity, gaussianBaselineConcentration, stochasticBaseline, stochasticStepSize, stochasticRate, plottingHowMuchFoodIsVisible, plottingWidth, plottingHeight, foodEstimatesPerUnitWorldLength, stochasticFoodScaling, responseScalingFactor, idealisedBugStepSize, realisticBugConstant, realisticBugLinear, realisticBugQuadratic, realisticBugMemoryLength, memoryLengthFactor, minimumAllowedMemoryLength, scaleFitnessBy, minimumPossibleDecimalValue, initialSensitivity, finalSensitivity, sensitivity, minimumConstantAllowed, maximumConstantAllowed, minimumLinearAllowed, maximumLinearAllowed, epsilon;
    public static boolean graphicalOutput, orkunTransition, printGeneration, printPreAndPostDurations, printBothPreAndPost, printCPUDurations, allowChangeOfTopology, differentBugsDifferentDrifts, useStochasticFood, stochasticIsFoodSquared, useIdealisedBug, equilibrateStochasticFood, calculateLinearFromQuadratic, useRealisticBug, calculateBugStepSize, loadNetworkPost, useConstantAndLinearForAdaptive, pretend2D, constrainSensitivity, constrainSensitivityInvertedEncodedByA, useNewWayOfMutatingAAndB, curveFitting, showCounter, showPlotLabel, showCorrelationTimeLength, showBug, instantaneousTumbles, averageStochasticFood;
    public static boolean[] responseFunction, fitnessFunction, proteinRoles, fitnessCalculation, stimulusCourse;
    public static List<Double> thresholds;
    public static String plotLabel;
    
    public Parameters(int jobNumber, int parameterSet) {
        mutationRateOriginal = new double[15];
        mutationRate = new double[mutationRateOriginal.length];
        responseFunction = new boolean[2];
        fitnessFunction = new boolean[5];
        fitnessCalculation = new boolean[4];
        stimulusCourse = new boolean[8];
        thresholds = new ArrayList<Double>();

        for(int i = 0; i < mutationRate.length; i++) {
            thresholds.add(i, 0.0);
        }

        Properties properties = new Properties();
        String filename = "parameters" + parameterSet;

        if(Modularity.whatTodo == 1) {
            filename = "./" + jobNumber + "/" + filename;
        }

		try {
			FileInputStream inputStream = new FileInputStream(filename);
			properties.load(inputStream);
		} catch(IOException e) {
            System.out.println("parameters file not present.");
			System.exit(-47);
		}

        //used in Modularity
        //1: functional form, 2: single_bug
        whatToSimulate = Integer.parseInt(properties.getProperty("whatToSimulate"));
        numberOfGenerations = Integer.parseInt(properties.getProperty("numberOfGenerations"));
        graphicalOutput = Boolean.parseBoolean(properties.getProperty("graphicalOutput"));
        printEvery = Integer.parseInt(properties.getProperty("printEvery"));
        numberOfConcurrentBugs = Integer.parseInt(properties.getProperty("numberOfConcurrentBugs"));
        howManyBugsToConsider = Integer.parseInt(properties.getProperty("howManyBugsToConsider"));
        maxDurationFactor = Double.parseDouble(properties.getProperty("maxDurationFactor"));
        maxWait = Integer.parseInt(properties.getProperty("maxWait"));
        foodRunningSumLimit = Integer.parseInt(properties.getProperty("foodRunningSumLimit"));
        flushEvery = Integer.parseInt(properties.getProperty("flushEvery"));
        printGeneration = Boolean.parseBoolean(properties.getProperty("printGeneration"));
        printGenerationEvery = Integer.parseInt(properties.getProperty("printGenerationEvery"));
        printPreAndPostDurations = Boolean.parseBoolean(properties.getProperty("printPreAndPostDurations"));
        loadNetwork = Integer.parseInt(properties.getProperty("loadNetwork"));
        loadNetworkPost = Boolean.parseBoolean(properties.getProperty("loadNetworkPost"));
        printBothPreAndPost = Boolean.parseBoolean(properties.getProperty("printBothPreAndPost"));
        printCPUDurations = Boolean.parseBoolean(properties.getProperty("printCPUDurations"));
        differentBugsDifferentDrifts = Boolean.parseBoolean(properties.getProperty("differentBugsDifferentDrifts"));
        useStochasticFood = Boolean.parseBoolean(properties.getProperty("useStochasticFood"));
        useIdealisedBug = Boolean.parseBoolean(properties.getProperty("useIdealisedBug"));
        useRealisticBug = Boolean.parseBoolean(properties.getProperty("useRealisticBug"));
        equilibrateStochasticFood = Boolean.parseBoolean(properties.getProperty("equilibrateStochasticFood"));
        idealisedBugStepSize = Double.parseDouble(properties.getProperty("idealisedBugStepSize"));
        curveFitting = Boolean.parseBoolean(properties.getProperty("curveFitting"));

        //used in Mutagenesis
        maxNumberOfProteins = Integer.parseInt(properties.getProperty("maxNumberOfProteins"));
        mutationRateOriginal[0] = Double.parseDouble(properties.getProperty("interactionMutationRate"));
        mutationRateOriginal[1] = Double.parseDouble(properties.getProperty("selfMutationRate"));
        mutationRateOriginal[2] = Double.parseDouble(properties.getProperty("sensitivityMutationRate"));
        mutationRateOriginal[3] = Double.parseDouble(properties.getProperty("proteinDuplicationRate"));
        mutationRateOriginal[4] = Double.parseDouble(properties.getProperty("proteinDeletionRate"));
        mutationRateOriginal[5] = Double.parseDouble(properties.getProperty("roleSwitchRate"));
        mutationRateOriginal[6] = Double.parseDouble(properties.getProperty("recruitmentRate"));
        mutationRateOriginal[7] = Double.parseDouble(properties.getProperty("effectorMotorAssociationRateRate"));
        mutationRateOriginal[8] = Double.parseDouble(properties.getProperty("effectorMotorDissociationRateRate"));
        mutationRateOriginal[9] = Double.parseDouble(properties.getProperty("responseScalingFactorRate"));
        mutationRateOriginal[10] = Double.parseDouble(properties.getProperty("realisticBugConstantRate"));
        mutationRateOriginal[11] = Double.parseDouble(properties.getProperty("realisticBugLinearRate"));
        mutationRateOriginal[12] = Double.parseDouble(properties.getProperty("realisticBugQuadraticRate"));
        mutationRateOriginal[13] = Double.parseDouble(properties.getProperty("realisticBugMemoryLengthRate"));
        mutationRateOriginal[14] = Double.parseDouble(properties.getProperty("epsilonRate"));
        deletionToAdditionRatio = Integer.parseInt(properties.getProperty("deletionToAdditionRatio"));
        maxInteractionStrength = Double.parseDouble(properties.getProperty("maxInteractionStrength"));
        minNumberOfProteins = Integer.parseInt(properties.getProperty("minNumberOfProteins"));
        allowChangeOfTopology = Boolean.parseBoolean(properties.getProperty("allowChangeOfTopology"));

        //used in ProteinNetwork
        initialProteinCount = Integer.parseInt(properties.getProperty("initialProteinCount"));
        responseScalingFactor = Double.parseDouble(properties.getProperty("responseScalingFactor"));

        //used in RungeKutta
        maxNumberOfIterations = Integer.parseInt(properties.getProperty("maxNumberOfIterations"));
        advanceTimeByDefault = Double.parseDouble(properties.getProperty("advanceTimeByDefault"));
        advanceTimeByMax = Double.parseDouble(properties.getProperty("advanceTimeByMax"));
        absoluteTolerance = Double.parseDouble(properties.getProperty("absoluteTolerance"));
        relativeTolerance = Double.parseDouble(properties.getProperty("relativeTolerance"));

        //used in World
        populationSize = Integer.parseInt(properties.getProperty("populationSize"));
        timeConstant = Double.parseDouble(properties.getProperty("timeConstant"));
        simulatedAnnealingTemperature = Double.parseDouble(properties.getProperty("simulatedAnnealingTemperature"));
        responseFunction[0] = Boolean.parseBoolean(properties.getProperty("evolveLinearResponse"));
        responseFunction[1] = Boolean.parseBoolean(properties.getProperty("evolveAdaptiveResponse"));
        fitnessFunction[0] = Boolean.parseBoolean(properties.getProperty("useHillClimbing"));
        fitnessFunction[1] = Boolean.parseBoolean(properties.getProperty("useKimura"));
        fitnessFunction[2] = Boolean.parseBoolean(properties.getProperty("useScaledKimura"));
        fitnessFunction[3] = Boolean.parseBoolean(properties.getProperty("useMetropolisHastings"));
        fitnessFunction[4] = Boolean.parseBoolean(properties.getProperty("useMetropolisHastingsStrict"));
        gaussianFunctionParameterC = Double.parseDouble(properties.getProperty("gaussianFunctionParameterC"));
        fitnessCalculation[0] = Boolean.parseBoolean(properties.getProperty("useLeastSquares"));
        fitnessCalculation[1] = Boolean.parseBoolean(properties.getProperty("useRewardingForFit"));
        fitnessCalculation[2] = Boolean.parseBoolean(properties.getProperty("useDecomposition"));
        fitnessCalculation[3] = Boolean.parseBoolean(properties.getProperty("useOrkun"));
        orkunAdaptiveWeighting = Double.parseDouble(properties.getProperty("orkunAdaptiveWeighting"));
        orkunLinearWeighting = Double.parseDouble(properties.getProperty("orkunLinearWeighting"));
        orkunTransition = Boolean.parseBoolean(properties.getProperty("orkunTransition"));
        orkunTransitionAt = Integer.parseInt(properties.getProperty("orkunTransitionAt"));
        orkunTransitionDuration = Integer.parseInt(properties.getProperty("orkunTransitionDuration"));
        stimulusCourse[0] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse1"));
        stimulusCourse[1] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse2"));
        stimulusCourse[2] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse3"));
        stimulusCourse[3] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse4"));
        stimulusCourse[4] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse5"));
        stimulusCourse[5] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse6"));
        stimulusCourse[6] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse7"));
        stimulusCourse[7] = Boolean.parseBoolean(properties.getProperty("useStimulusCourse8"));
        baselineFitness = Double.parseDouble(properties.getProperty("baselineFitness"));
        steadyStateThreshold = Double.parseDouble(properties.getProperty("steadyStateThreshold"));
        steadyStateQueueSize = Integer.parseInt(properties.getProperty("steadyStateQueueSize"));
        steadyStateMaxEstimates = Integer.parseInt(properties.getProperty("steadyStateMaxEstimates"));
        worldLength = Double.parseDouble(properties.getProperty("worldLength"));
        stepSize = Double.parseDouble(properties.getProperty("stepSize"));
        stochasticBaseline = Double.parseDouble(properties.getProperty("stochasticBaseline"));
        stochasticStepSize = Double.parseDouble(properties.getProperty("stochasticStepSize"));
        numberOfModes = Integer.parseInt(properties.getProperty("numberOfModes"));
        stochasticRate = Double.parseDouble(properties.getProperty("stochasticRate"));
        stochasticFoodScaling = Double.parseDouble(properties.getProperty("stochasticFoodScaling"));
        stochasticIsFoodSquared = Boolean.parseBoolean(properties.getProperty("stochasticIsFoodSquared"));
        scaleFitnessBy = Double.parseDouble(properties.getProperty("scaleFitnessBy"));
        averageStochasticFood = Boolean.parseBoolean(properties.getProperty("averageStochasticFood"));

        //used in Bug
        effectorMotorAssociationRate = Double.parseDouble(properties.getProperty("effectorMotorAssociationRate"));
        effectorMotorDissociationRate = Double.parseDouble(properties.getProperty("effectorMotorDissociationRate"));
        driftVelocity = Double.parseDouble(properties.getProperty("driftVelocity"));
        realisticBugConstant = Double.parseDouble(properties.getProperty("realisticBugConstant"));
        realisticBugLinear = Double.parseDouble(properties.getProperty("realisticBugLinear"));
        realisticBugQuadratic = Double.parseDouble(properties.getProperty("realisticBugQuadratic"));
        realisticBugMemoryLength = Double.parseDouble(properties.getProperty("realisticBugMemoryLength"));
        calculateLinearFromQuadratic = Boolean.parseBoolean(properties.getProperty("calculateLinearFromQuadratic"));
        memoryLengthFactor = Double.parseDouble(properties.getProperty("memoryLengthFactor"));
        calculateBugStepSize = Boolean.parseBoolean(properties.getProperty("calculateBugStepSize"));
        minimumAllowedMemoryLength = Double.parseDouble(properties.getProperty("minimumAllowedMemoryLength"));
        useConstantAndLinearForAdaptive = Boolean.parseBoolean(properties.getProperty("useConstantAndLinearForAdaptive"));
        pretend2D = Boolean.parseBoolean(properties.getProperty("pretend2D"));
        minimumPossibleDecimalValue = Double.parseDouble(properties.getProperty("minimumPossibleDecimalValue"));
        constrainSensitivity = Boolean.parseBoolean(properties.getProperty("constrainSensitivity"));
        constrainSensitivityInvertedEncodedByA = Boolean.parseBoolean(properties.getProperty("constrainSensitivityInvertedEncodedByA"));
        initialSensitivity = Double.parseDouble(properties.getProperty("initialSensitivity"));
        finalSensitivity = Double.parseDouble(properties.getProperty("finalSensitivity"));
        startIncreasingSensitivityAtGeneration = Integer.parseInt(properties.getProperty("startIncreasingSensitivityAtGeneration"));
        useNewWayOfMutatingAAndB = Boolean.parseBoolean(properties.getProperty("useNewWayOfMutatingAAndB"));
        minimumConstantAllowed = Double.parseDouble(properties.getProperty("minimumConstantAllowed"));
        maximumConstantAllowed = Double.parseDouble(properties.getProperty("maximumConstantAllowed"));
        minimumLinearAllowed = Double.parseDouble(properties.getProperty("minimumLinearAllowed"));
        maximumLinearAllowed = Double.parseDouble(properties.getProperty("maximumLinearAllowed"));
        epsilon = Double.parseDouble(properties.getProperty("epsilon"));
        instantaneousTumbles = Boolean.parseBoolean(properties.getProperty("instantaneousTumbles"));

        if(constrainSensitivity || constrainSensitivityInvertedEncodedByA) {
            sensitivity = initialSensitivity;
        } else if(mutationRateOriginal[14] > 0) {
            sensitivity = Double.parseDouble(properties.getProperty("sensitivity"));
        } else {
            sensitivity = Double.MAX_VALUE;
        }

        if(stimulusCourse[0]) {
            stimulusCourseDuration = 150;
        } else if(stimulusCourse[1]) {
            stimulusCourseDuration = 240;
        } else if(stimulusCourse[2]) {
            stimulusCourseDuration = 3.54;
        } else if(stimulusCourse[3]) {
            stimulusCourseDuration = 875;
        } else if(stimulusCourse[4]) {
            stimulusCourseDuration = 3000;
        } else if(stimulusCourse[5]) {
            //this is just a dummy number - it is filled in ProteinNetwork/loadRealisticBug
            stimulusCourseDuration = 1100;
        } else if(stimulusCourse[6]) {
            stimulusCourseDuration = 2 / stochasticRate;
        } else if(stimulusCourse[7]) {
            stimulusCourseDuration = 1625;
        }

        //used in Protein
        int j = 0;

        for(int i = 1; i <= 100; i++) {
            if(properties.getProperty("isProtein" + i + "AKinase") != null) {
                j++;
            }
            else {
                break;
            }
        }

        proteinRoles = new boolean[j];

        for(int i = 1; i <= j; i++) {
            proteinRoles[i-1] = Boolean.parseBoolean(properties.getProperty("isProtein" + i + "AKinase"));
        }

        //used in Gaussian
        gaussianDetectionThreshold = Double.parseDouble(properties.getProperty("gaussianDetectionThreshold"));
        gaussianBaselineConcentration = Double.parseDouble(properties.getProperty("gaussianBaselineConcentration"));
        j = 0;

        for(int i = 1; i <= 100; i++) {
            if(properties.getProperty("gaussian" + i + "Width") != null) {
                j++;
            } else {
                break;
            }
        }

        gaussianWidths = new double[j];
        gaussianHeights = new double[j];
        gaussianCentres = new double[j];

        for(int i = 0; i < j; i++) {
            gaussianWidths[i] = Double.parseDouble(properties.getProperty("gaussian" + (i+1) + "Width"));
            gaussianHeights[i] = Double.parseDouble(properties.getProperty("gaussian" + (i+1) + "Height"));
            gaussianCentres[i] = Double.parseDouble(properties.getProperty("gaussian" + (i+1) + "Centre"));
        }

        //used in PlottingStochastic and/or PlottingSingleBug
        plottingHowMuchFoodIsVisible = Double.parseDouble(properties.getProperty("plottingHowMuchFoodIsVisible"));
        plottingWidth = Double.parseDouble(properties.getProperty("plottingWidth"));
        plottingHeight = Double.parseDouble(properties.getProperty("plottingHeight"));
        foodEstimatesPerUnitWorldLength = Double.parseDouble(properties.getProperty("foodEstimatesPerUnitWorldLength"));
        showCounter = Boolean.parseBoolean(properties.getProperty("showCounter"));
        showPlotLabel = Boolean.parseBoolean(properties.getProperty("showPlotLabel"));
        plotLabel = properties.getProperty("plotLabel");
        showCorrelationTimeLength = Boolean.parseBoolean(properties.getProperty("showCorrelationTimeLength"));
        showBug = Boolean.parseBoolean(properties.getProperty("showBug"));
    }

    public void checkConsistency() {
        if((responseFunction[0] && responseFunction[1]) || (!responseFunction[0] && !responseFunction[1])) {
            System.out.println("Inconsistency in the parameters file: unclear which response function the network should evolve to imitate.");
            System.exit(-47);
        }

        int fitnessFunctionCount = 0, fitnessCalculationCount = 0, stimulusCourseCount = 0;

        for(Boolean b : fitnessFunction) {
            if(b) {
                fitnessFunctionCount++;
            }
        }
        if(fitnessFunctionCount != 1) {
            System.out.println("Inconsistency in the parameters file: unclear which fitness function should be used.");
            System.exit(-47);
        }

        for(Boolean b : fitnessCalculation) {
            if(b) {
                fitnessCalculationCount++;
            }
        }
        if(fitnessCalculationCount != 1) {
            System.out.println("Inconsistency in the parameters file: unclear how fitness should be calculated.");
            System.exit(-47);
        }

        for(Boolean b : stimulusCourse) {
            if(b) {
                stimulusCourseCount++;
            }
        }
        if(stimulusCourseCount != 1) {
            System.out.println("Inconsistency in the parameters file: unclear which stimulus course should be used.");
            System.exit(-47);
        }

        if(initialProteinCount != proteinRoles.length) {
            System.out.println("Inconsistency in the parameters file: unclear how many proteins the initial network should contain.");
            System.exit(-47);
        }

        if(fitnessFunction[1] || fitnessFunction[2]) {
            System.out.println("Kimura and scaled Kimura won't work properly because of the way fitness is defined (when the fitness of the pre-mutagenesis network is negative, networks that fit the goal less well are selected for).");
            System.exit(-47);
        }

        if(fitnessCalculation[3]) {
            if(orkunLinearWeighting == 0 && orkunAdaptiveWeighting == 0) {
                System.out.println("Inconsistency in the parameters file: Orkun's fitness fitness should be used but both orkunLinearWeighting and orkunAdaptiveWeighting are set to 0.");
                System.exit(-47);
            }
        }

        if(orkunTransition) {
            if(orkunLinearWeighting == 0 || orkunAdaptiveWeighting == 0) {
                System.out.println("Inconsistency in the parameters file: orkunTransition is set to true but orkunLinearWeighting or orkunAdaptiveWeighting is set to 0.");
                System.exit(-47);
            } if(orkunTransitionAt + orkunTransitionDuration > numberOfGenerations) {
                System.out.println("Inconsistency in the parameters file: the sum of orkunTransitionAt and orkunTransitionDuration is larger than numberOfGenerations.");
                System.exit(-47);
            }
        }

        if(Modularity.whatTodo == 9 && whatToSimulate == 1) {
            System.out.println("Inconsistency: whatToDo = 9 but whatToSimulate = 1 (if whatToDo = 9, whatToSimulate should = 2).");
            System.exit(-47);
        }
    }

    public static void normaliseMutationRates(int proteinCount, int receptorCount) {
        double sum = 0;

        mutationRate[0] = proteinCount * (proteinCount - 1) * mutationRateOriginal[0];
        mutationRate[1] = 2 * proteinCount * mutationRateOriginal[1]; //every protein has a passive activation and deactivation rate which is why the proteinCount is multiplied by 2 (2 rates)
        /* three things we might wanna do with sensitivities:
           1) not mutate them and vary sensitivity manually - in that case it doesn't matter how mutationRate[2] is defined because mutationRateOriginal[2] = 0
           2) when a mutation occurs, the change is applied to all the sensitivities (associated with all the receptors) - not biologically relevant since in other aspects we consider the receptors to be genetically distinct
           3) mutate them "independently" (as opposed to 2) )
           if we choose to implement 2), we should not multiply mutationRateOriginal by receptorCount since sensitivity is, in reality, an entity separable from receptors so the number of receptors has nothing to do with it.
           if, however, we choose to implement 3), mutationRate SHOULD be defined the way it is so as to reflect that each of the sensitivities is a separate entity.
         */
        mutationRate[2] = receptorCount * mutationRateOriginal[2];
        mutationRate[3] = proteinCount * mutationRateOriginal[3];
        mutationRate[4] = proteinCount * mutationRateOriginal[4];
        mutationRate[5] = proteinCount * mutationRateOriginal[5];
        mutationRate[6] = 2 * proteinCount * mutationRateOriginal[6];
        mutationRate[7] = mutationRateOriginal[7];
        mutationRate[8] = mutationRateOriginal[8];
        mutationRate[9] = mutationRateOriginal[9];
        mutationRate[10] = mutationRateOriginal[10];
        mutationRate[11] = mutationRateOriginal[11];
        mutationRate[12] = mutationRateOriginal[12];
        mutationRate[13] = mutationRateOriginal[13];
        mutationRate[14] = mutationRateOriginal[14];

        for(double d : mutationRate)
            sum += d;

        for(int i = 0; i < mutationRate.length; i++)
            mutationRate[i] /= sum;

        for(int i = 0; i < mutationRate.length; i++) {
            if(i == 0) thresholds.set(0, mutationRate[0]);
            else thresholds.set(i, thresholds.get(i-1) + mutationRate[i]);
        }
    }
}
