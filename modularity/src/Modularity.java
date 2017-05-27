import com.google.common.base.Joiner;
import org.apache.commons.math3.random.BitsStreamGenerator;
import org.apache.commons.math3.random.Well44497b;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.*;
import javax.swing.*;

public class Modularity {
    public static Parameters parameters;
    public static PrintWriter output = null, fitnesses = null, orkunComponents = null, network = null, interactions = null, receptors = null, roles = null, effectors = null, fitnessRejectedWriter = null, fixationProbabilities = null, effectorMotorRates = null, durations = null, positions = null, tumblingProbabilities = null, runDurations = null, tumbleDurations = null, alphas = null, runningProbabilities = null, convolutions = null;
    public static DataInputStream urandom = null;
    private List<String> filenames;
    //1: cluster, 2: master, 3: output_analysis, 4: testRandomBug, 5: testRKProteinNetwork, 6: testRKExponentialDecay, 7: testMutagenesis, 8: test, 9: testSingleBugOneGeneration, 10: testStochasticFood
    public static int whatTodo;
    private int loadGeneration;
    public static BitsStreamGenerator random;
    public static List<List<Double>> packer;

    public Modularity(int jobNumber, int parameterSet, int whatTodo, int run, int loadGeneration) {
        Modularity.whatTodo = whatTodo;
        parameters = new Parameters(jobNumber, parameterSet);
        parameters.checkConsistency();
        filenames = new ArrayList<String>();
        this.loadGeneration = loadGeneration;

        if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2 || Modularity.whatTodo == 3 || Modularity.whatTodo == 9) {
            initialiseFilenames(jobNumber, parameterSet, run);
            clearFiles();
            initialiseFiles();
        }

        long seed = getSeed();
        random = new Well44497b(seed);
        System.out.format("%d\n", seed);
        packer = new ArrayList<List<Double>>();
        for(int i = 0; i < 2; i++) packer.add(new ArrayList<Double>());
    }

    public static void main(String[] args) {
        Modularity modularity;

        //0: number of parameters files, 1: how many runs per each parameters file, 2: output files from which run to load (output analysis)/which folder to load the parameters file from and output the files to (cluster), 3: what to do, 4: load from the effector_motor_rates file the effector-motor binding rates from the specified generation (if different from -1)
        for(int parameterSet = 1; parameterSet <= Integer.parseInt(args[0]); parameterSet++) {
            System.out.println("Parameter set " + parameterSet + " (of " + Integer.parseInt(args[0]) + ")");

            for(int run = 1; run <= Integer.parseInt(args[1]); run++) {
                System.out.println("Run " + run + " (of " + Integer.parseInt(args[1]) + ")");
                modularity = new Modularity(Integer.parseInt(args[2]), parameterSet, Integer.parseInt(args[3]), run, Integer.parseInt(args[4]));

                if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2) {
                    if(Parameters.whatToSimulate == 1) modularity.runSimulation();
                    if(Parameters.whatToSimulate == 2) modularity.singleBugEvolution();
                } else if(Modularity.whatTodo == 3) {
                    modularity.seeOutput();
                } else if(Modularity.whatTodo == 4) {
                    modularity.testRandomBug();
                } else if(Modularity.whatTodo == 5) {
                    modularity.testRKProteinNetwork();
                } else if(Modularity.whatTodo == 6) {
                    modularity.testRKExponentialDecay();
                } else if(Modularity.whatTodo == 7) {
                    modularity.testMutagenesis();
                } else if(Modularity.whatTodo == 8) {
                    modularity.test();
                } else if(Modularity.whatTodo == 9) {
                    modularity.testSingleBugOneGeneration();
                } else if(Modularity.whatTodo == 10) {
                    modularity.testStochasticFood();
                }

                if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2 || Modularity.whatTodo == 3 || Modularity.whatTodo == 9) {
                    modularity.closeFiles();
                }
            }
        }

        System.out.println("Done.");
        System.exit(0);
    }

    private void runSimulation() {
        Plotting plot = new Plotting();
        int acceptanceProbability = -1;
        ProteinNetwork proteinNetworkPre = new ProteinNetwork();
        proteinNetworkPre.initialise();
        RungeKutta rungeKuttaPre = new RungeKutta(null);
        RungeKutta rungeKuttaPost = null;
        World worldPre = new World(proteinNetworkPre.getProteinIds());
        World worldPost;
        World.setSensitivity();
        proteinNetworkPre.setSensitivity(World.getSensitivity());
        ProteinNetwork proteinNetworkPost = null;
        Mutagenesis mutagenesis = new Mutagenesis();
        mutagenesis.fillNetwork(proteinNetworkPre);
        double fitnessRejected = -47;
        Bug bug = new Bug(-1);

        if(Modularity.whatTodo == 2 && Parameters.graphicalOutput) {
            plot.initialise();
        }

        for(int generation = 1; generation <= Parameters.numberOfGenerations; generation++) {
            if(Parameters.printGeneration && generation % Parameters.printGenerationEvery == 0) {
                System.out.println("Generation " + generation);
            }

            if(Parameters.whatToSimulate == 1 && Parameters.fitnessCalculation[3] && Parameters.orkunTransition) {
                World.transitionFromLinearToAdaptive(generation);
            }

            try { proteinNetworkPost = (ProteinNetwork) proteinNetworkPre.clone();
            } catch(CloneNotSupportedException e) { System.exit(-47); }
            try { rungeKuttaPost = (RungeKutta) rungeKuttaPre.clone();
            } catch(CloneNotSupportedException e) { System.exit(-47); }

            mutagenesis.mutate(proteinNetworkPost, rungeKuttaPost);
            //proteinNetworkPre.print();
            //proteinNetworkPost.print();

            if(generation == 1) {
                worldPre = new World(proteinNetworkPre.getProteinIds());
                worldPre.waitUntilSteadyState(rungeKuttaPre, proteinNetworkPre, bug);

                if(worldPre.isEquilibrated()) {
                    while(worldPre.getTime() <= Parameters.stimulusCourseDuration) {
                        if(!rungeKuttaPre.cashKarp(proteinNetworkPre, worldPre, bug, false)) {
                            System.out.println("fourthOrder used!");
                            System.exit(-47);
                            //rungeKuttaPre.fourthOrder(proteinNetworkPre, worldPre.getInput());
                        }

                        worldPre.increaseTime(rungeKuttaPre.getStepSizeThisStep());
                        worldPre.logTime();
                        worldPre.addEffectorActivity(rungeKuttaPre.getActivityLevel(Protein.EXTRA_INTERACTIONS + 1));
                        worldPre.advanceInput();

                        if(!Parameters.fitnessCalculation[3]) {
                            worldPre.calculateGoal();
                        }
                    }

                    worldPre.dropEstimates(null);
                    //worldPre.calculateAverageDesired();
                    worldPre.calculateFitness();
                    //System.out.println("goodnessOfFitPre: " + worldPre.getGoodnessOfFit());
                } else {
                    worldPre.setFitness(0);
                }

                rungeKuttaPre.resetActivityLevels();
            }

            if(Parameters.fitnessCalculation[3] && Parameters.orkunTransition && generation >= Parameters.orkunTransitionAt && generation <= Parameters.orkunTransitionAt + Parameters.orkunTransitionDuration) {
                worldPre.calculateFitness();
            }

            worldPost = new World(proteinNetworkPost.getProteinIds());
            worldPost.waitUntilSteadyState(rungeKuttaPost, proteinNetworkPost, bug);

            if(worldPost.isEquilibrated()) {
                while(worldPost.getTime() <= Parameters.stimulusCourseDuration) {
                    if(!rungeKuttaPost.cashKarp(proteinNetworkPost, worldPost, bug, false)) {
                        System.out.println("fourthOrder used!");
                        System.exit(-47);
                        //rungeKuttaPost.fourthOrder(proteinNetworkPost, worldPost.getInput());
                    }

                    worldPost.increaseTime(rungeKuttaPost.getStepSizeThisStep());
                    worldPost.logTime();
                    worldPost.addEffectorActivity(rungeKuttaPost.getActivityLevel(Protein.EXTRA_INTERACTIONS + 1));
                    worldPost.advanceInput();

                    if(!Parameters.fitnessCalculation[3]) {
                        worldPost.calculateGoal();
                    }
                }

                worldPost.dropEstimates(null);
                //worldPost.calculateAverageDesired();
                worldPost.calculateFitness();
                //System.out.println("goodnessOfFitPost: " + worldPost.getGoodnessOfFit());
            } else {
                worldPost.setFitness(0);
            }

            rungeKuttaPost.resetActivityLevels();

            if(Parameters.fitnessCalculation[2]) {
                if(Parameters.responseFunction[0]) {
                    acceptanceProbability = World.fitnessFunction(worldPre.getFitness(0), worldPost.getFitness(0));
                } else if(Parameters.responseFunction[1]) {
                    acceptanceProbability = World.fitnessFunction(worldPre.getFitness(1), worldPost.getFitness(1));
                }
            } else {
                acceptanceProbability = World.fitnessFunction(worldPre.getFitness(World.DEFAULT_FITNESS), worldPost.getFitness(World.DEFAULT_FITNESS));
            }

            if(acceptanceProbability == 1) {
                fitnessRejected = worldPre.getFitness();

                try { proteinNetworkPre = (ProteinNetwork) proteinNetworkPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }
                try { rungeKuttaPre = (RungeKutta) rungeKuttaPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }
                try { worldPre = (World) worldPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }

                worldPre.setOutputsSteadyState(worldPost.cloneOutputsSteadyState());
            } else if(acceptanceProbability == 0) {
                fitnessRejected = worldPost.getFitness();
            } else if(acceptanceProbability == -1) {
                System.out.println("An error has been encountered while executing the fitness function.");
                System.exit(-47);
            }

            //System.out.println("linear: " + worldPre.getFitness(0) + ", adaptive: " + worldPre.getFitness(1) + ", flat: " + worldPre.getFitness(2));

            if(networkShouldBePrinted(generation)) proteinNetworkPre.appendToFile(generation);
            worldPre.appendToFile(generation, networkShouldBePrinted(generation), fitnessRejected);

            if(Parameters.flushEvery != -1 && generation % Parameters.flushEvery == 0) {
                flush();
            }

            if(Parameters.fitnessFunction[2]) {
                World.setGamma(generation);
            }

            rungeKuttaPre.resetStepSizeNextStep();

            if(Modularity.whatTodo == 2 && Parameters.graphicalOutput && generation % 10 == 0 && worldPre.isEquilibrated()) {
                plot.clearTraces();

                //It's important that the points are added in this order i.e. steady state first because otherwise the list looks like this: (just the times) 1, 2, 3, ..., 149, 150, -2, -1, 0 and a line is drawn between the point 150,y1 and -2,y2.
                for(int i = 0; i < worldPre.getSteadyStateDataPointCount(); i++) {
                    plot.updateTrace(0, worldPre.getTimesSteadyState(i), worldPre.getOutputsSteadyState(Protein.EXTRA_INTERACTIONS + 1, i));
                    plot.updateTrace(1, worldPre.getTimesSteadyState(i), worldPre.getInputsSteadyState(i));
                }
                for(int i = 0; i < worldPre.getDataPointCount(); i++) {
                    plot.updateTrace(0, worldPre.getTime(i), worldPre.getEffectorActivity(i));
                    plot.updateTrace(1, worldPre.getTime(i), worldPre.getInput(i));

                    if(!Parameters.fitnessCalculation[3]) {
                        plot.updateTrace(3, worldPre.getTime(i), worldPre.getGoal(i));
                    }
                }
            }

            //System.out.println("----------------------------------------------------------------------------------------------");
        }

        if(Modularity.whatTodo == 2 && Parameters.graphicalOutput) {
            plot.destroy();
        }
    }

    private void testRandomBug() {
        Bug bug = new Bug(-1);
        JFrame jFrame = new JFrame("Chemotaxis");
        World.initialiseGaussians();
        PlottingSingleBug plot = new PlottingSingleBug(bug.getPosition(), World.getGaussians());

        if(Parameters.graphicalOutput) {
            jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
            jFrame.add(plot);
            jFrame.pack();
            jFrame.setVisible(true);
        }

        for(int i = 0; i < 10000; i++) {
            bug.nowWhatRandom();

            if(Parameters.graphicalOutput && !bug.isTumbling()) {
                plot.moveBug(bug.getPosition());
            }

            try {
                Thread.sleep(10);
            } catch(InterruptedException e) {}
        }

        jFrame.dispose();
    }

    public boolean networkShouldBePrinted(int generation) {
        if(generation % Parameters.printEvery == 0) {
            return true;
        } else if(Parameters.orkunTransition) {
            if(generation == Parameters.orkunTransitionAt - 1) {
                return true;
            } else if(generation > Parameters.orkunTransitionAt && generation < Parameters.orkunTransitionAt + 40000 && generation % 50 == 0) {
                return true;
            } else if(generation > Parameters.orkunTransitionAt + 40000 && generation % 100 == 0) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    private void seeOutput() {
        NumberFormat format = new DecimalFormat("0.######E0");
        ProteinNetwork proteinNetwork = loadNetwork();
        proteinNetwork.loadRates();

        if(Parameters.useIdealisedBug) {
            proteinNetwork.loadResponseScalingFactor();
        } if(Parameters.useRealisticBug) {
            proteinNetwork.loadRealisticBug();
        }

        if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            proteinNetwork.fillTauLogSpace();
        }

        RungeKutta rk = new RungeKutta(proteinNetwork.getProteinIds());
        World world = new World(proteinNetwork.getProteinIds());
        world.setRandom(Modularity.random);
        world.memoryCalculations(proteinNetwork);

        if(Parameters.useStochasticFood && Parameters.equilibrateStochasticFood) {
            World.determineEquilibrationPeriod(proteinNetwork.realisticBugParameters.get(3), proteinNetwork.realisticBugParameters.get(3));
            world.equilibrateStochasticFood();
        }

        Map<Integer, List<Double>> outputs = new HashMap<Integer, List<Double>>();
        double food = 0, limit = 0, currentFood = 0, startedRunning = 0, startedTumbling = 0, probTumbling, foodMidStep, stepSize = Parameters.useRealisticBug ? proteinNetwork.bugStepSize : Parameters.idealisedBugStepSize;
        int index = -1, indexOld = -1;
        long timeBefore = System.nanoTime();
        boolean isTumblingPast = false;

        for(Integer i : proteinNetwork.getProteinIds()) {
            outputs.put(i, new ArrayList<Double>());
        }

        Bug bug = new Bug(Parameters.driftVelocity);
        bug.setRandom(Modularity.random);

        if(Parameters.useIdealisedBug) {
            bug.initialiseFoodPast(world.stochasticFoodDistribution(bug.getPosition(), -1));
        } if(Parameters.useRealisticBug) {
            bug.initialiseFoodsPast(world.getFoodsPast(proteinNetwork, bug));
            world.calculateExponentials(proteinNetwork);
            bug.setExponentials(world.getExponentials());
            bug.setLocalCopiesOfParameters(proteinNetwork);
            bug.equilibrateAlpha();
        }

        World.initialiseGaussians();
        List<Double> foodDistribution = new ArrayList<Double>();

        for(double d = 0; d <= Parameters.worldLength; d += 1 / Parameters.foodEstimatesPerUnitWorldLength) {
            foodDistribution.add(world.stochasticFoodDistribution(d, -1));
        }

        SaveImages saveImages = new SaveImages(bug.getPosition(), foodDistribution, world.getTime());
        PlottingSingleBug plot = new PlottingSingleBug(bug.getPosition(), World.getGaussians());
        PlottingStochastic plotStochastic = new PlottingStochastic(bug.getPosition(), foodDistribution);
        JFrame jFrame = new JFrame("Chemotaxis");

        if(Parameters.whatToSimulate == 2 && Parameters.graphicalOutput) {
            jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

            if(Parameters.useStochasticFood) {
                jFrame.add(plotStochastic);
            } else {
                jFrame.add(plot);
            }

            jFrame.pack();
            jFrame.setVisible(true);
        }

        if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
            world.waitUntilSteadyState(rk, proteinNetwork, bug);
        }

        while(itIsNotTimeToStopTheRun(world, rk, proteinNetwork)) {
            world.increaseSteps();

            if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
                if(!rk.cashKarp(proteinNetwork, world, bug, false)) {
                    System.out.println("fourthOrder used!");
                    System.exit(-47);
                }

                world.increaseTime(rk.getStepSizeThisStep());
            } else {
                world.setTime(world.getSteps() * stepSize);
            }

            if(Parameters.whatToSimulate == 1) {
                world.logTime();

                if(Parameters.useIdealisedBug) {
                    world.addEffectorActivity(bug.moveIdealised(proteinNetwork, Parameters.idealisedBugStepSize, World.getInput(world.getTime())));
                } else if(Parameters.useRealisticBug) {
                    world.addEffectorActivity(bug.moveRealisticFasterTetheredBug(World.getInput(world.getTime() - stepSize), World.getInput(world.getTime())));
                } else {
                    world.addEffectorActivity(rk.getActivityLevel(Protein.EXTRA_INTERACTIONS + 1));

                    for(Integer i : proteinNetwork.getProteinIds()) {
                        outputs.get(i).add(rk.getActivityLevel(i));
                    }
                }

                world.advanceInput();
                world.calculateGoal();
            } else if(Parameters.whatToSimulate == 2) {
                if(Parameters.graphicalOutput) {
                    if(Parameters.useStochasticFood) {
                        if(index != indexOld) {
                            foodDistribution.clear();

                            for(double d = 0; d <= Parameters.worldLength; d += 1 / Parameters.foodEstimatesPerUnitWorldLength) {
                                foodDistribution.add(world.stochasticFoodDistribution(d, index));
                            }

                            plotStochastic.moveFood(foodDistribution);
                            plotStochastic.moveFoodRepaint();
                            saveImages.moveFood(foodDistribution);
                        }

                        plotStochastic.moveBug(bug.getPosition());
                        plotStochastic.moveBugRepaint();

                        if(world.isTimeWhole(proteinNetwork.howManyStepsInUnitTime)) {
                            saveImages.changeTime(world.getTime());
                            saveImages.moveBug(bug.getPosition());
                            saveImages.saveImage();
                        }
                    } else {
                        plot.moveBug(bug.getPosition());
                    }
                }

                indexOld = index;
                index = World.getIndexOfNearestStochasticTime(world.getTime());
                foodMidStep = world.stochasticFoodDistribution(bug.getMidPosition(), World.getIndexOfNearestStochasticTime(world.getTime() - stepSize));
                currentFood = world.stochasticFoodDistribution(bug.getPosition(), index);

                if(Parameters.useIdealisedBug) {
                    bug.moveIdealised(proteinNetwork, Parameters.idealisedBugStepSize, currentFood);
                } else if(Parameters.useRealisticBug) {
                    if(!Parameters.instantaneousTumbles) {
                        probTumbling = bug.moveRealisticFaster(foodMidStep, currentFood);
                    } else {
                        probTumbling = bug.moveInstantaneousTumbles(foodMidStep, currentFood);
                    }

                    if(world.isTimeWhole(proteinNetwork.howManyStepsInUnitTime * 10)) {
                        if(!isTumblingPast) {
                            tumblingProbabilities.format("%f,%s\n", world.getTime(), format.format(probTumbling));
                        }
                        if(isTumblingPast) {
                            runningProbabilities.format("%f,%s\n", world.getTime(), format.format(probTumbling));
                        }

                        alphas.format("%f,%s\n", world.getTime(), format.format(bug.getAlpha()));
                        convolutions.format("%f,%s\n", world.getTime(), format.format(bug.getConvolution()));
                    }

                    if(!isTumblingPast && bug.isTumbling()) {
                        runDurations.format("%f\n", world.getTime() - startedRunning);
                        startedTumbling = world.getTime();
                        isTumblingPast = true;
                    }
                    if(isTumblingPast && !bug.isTumbling()) {
                        tumbleDurations.format("%f\n", world.getTime() - startedTumbling);
                        startedRunning = world.getTime();
                        isTumblingPast = false;
                    }
                } else {
                    bug.move(proteinNetwork, rk.getActivityLevel(Protein.EXTRA_INTERACTIONS + 1), rk.getStepSizeThisStep());
                }

                if(world.isTimeWhole(proteinNetwork.howManyStepsInUnitTime)) {
                    world.logTime();

                    if(Parameters.useStochasticFood) {
                        food += currentFood;
                        world.addInput(currentFood);
                    } else {
                        world.logFood(bug.getPosition());
                    }

                    if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
                        for(Integer i : proteinNetwork.getProteinIds()) {
                            outputs.get(i).add(rk.getActivityLevel(i));
                        }
                    }

                    if(world.getTime() % 100 == 0) {
                        System.out.println("Time: " + world.getTime() + ", food: " + food);
                    }

                    if(!Parameters.useStochasticFood) {
                        food += world.getFoodAt(bug.getPosition()) - Parameters.gaussianBaselineConcentration;
                    }

                    if(food > Parameters.foodRunningSumLimit) break;
                    positions.format("%f,%f\n", world.getTime(), bug.getPosition());
                    positions.flush();
                }
            }
        }

        System.out.println("CPU duration: " + (((double) (System.nanoTime() - timeBefore)) / 1000000000) + ", number of RK estimates: " + world.getSteps() + ", number of required RK estimates: " + world.getTime());
        jFrame.dispose();

        if(Parameters.whatToSimulate == 1) {
            if(!Parameters.useRealisticBug) {
                world.dropEstimates(outputs);
            }

            world.calculateFitness();
            System.out.println("Fitness: " + world.getFitness());
        }

        if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            output.format(",\"Protein4\"");
        } else {
            for(Integer i : proteinNetwork.getProteinIds()) {
                output.format(",\"Protein" + i + "\"");
            }
        }

        output.println();

        if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
            for(int i = 0; i < world.getSteadyStateDataPointCount(); i++) {
                output.format("%f,%f", world.getTimesSteadyState(i), world.getInputsSteadyState(i));

                for(Integer j : proteinNetwork.getProteinIds()) {
                    output.format(",%f", world.getOutputsSteadyState(j, i));
                }

                output.println();
            }
        }

        int upperLimit;

        if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            upperLimit = world.getDataPointCount();
        } else {
            upperLimit = outputs.get(Protein.EXTRA_INTERACTIONS + 1).size();
        }

        int step = 1;

        if(Parameters.stimulusCourse[4]) {
            step = 10;
        } else if(Parameters.stimulusCourse[5]) {
            step = 1;
        }

        for(int i = 0; i < upperLimit; i += step) {
            output.format("%f,%f", world.getTime(i), world.getInput(i));

            if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
                output.format(",%s", format.format(world.getEffectorActivity(i)));
            } else {
                for(Integer j : proteinNetwork.getProteinIds()) {
                    output.format(",%f", outputs.get(j).get(i));
                }
            }

            output.println();
        }

        if(Parameters.whatToSimulate == 1 && Parameters.fitnessCalculation[3] && !Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
            for(int i = 0; i < outputs.get(4).size() - 1; i++) {
                orkunComponents.format("%f,%f,%f,%f\n", world.getTime(i+1), world.getOrkunFitnessComponents(0, i), world.getOrkunFitnessComponents(1, i), world.getOrkunFitnessComponents(2, i));
            }
        }
    }

    private boolean itIsNotTimeToStopTheRun(World world, RungeKutta rk, ProteinNetwork proteinNetwork) {
        if(Parameters.whatToSimulate == 1) {
            return world.getTime() + rk.getStepSizeNextStep() <= Parameters.stimulusCourseDuration;
        } else if(Parameters.whatToSimulate == 2) {
            return world.getSteps() / proteinNetwork.howManyStepsInUnitTime < 10000;
        } else {
            return false;
        }
    }

    private void flush() {
        if(output != null) output.flush();
        if(fitnesses != null) fitnesses.flush();
        if(orkunComponents != null) orkunComponents.flush();
        if(network != null) network.flush();
        if(interactions != null) interactions.flush();
        if(receptors != null) receptors.flush();
        if(roles != null) roles.flush();
        if(effectors != null) effectors.flush();
        if(fitnessRejectedWriter != null) fitnessRejectedWriter.flush();
        if(fixationProbabilities != null) fixationProbabilities.flush();
        if(effectorMotorRates != null) effectorMotorRates.flush();
        if(durations != null) durations.flush();
    }

    private void singleBugEvolution() {
        int acceptanceProbability;
        double fitnessRejected = -47;
        ProteinNetwork proteinNetworkPre, proteinNetworkPost = null;
        long timeBefore, timeStart = System.nanoTime(), seed;

        if(Parameters.curveFitting) {
            loadCSV("packer.csv", packer);
        }

        if(Parameters.loadNetwork == -1) {
            proteinNetworkPre = new ProteinNetwork();
            proteinNetworkPre.initialise();
        } else {
            proteinNetworkPre = loadNetwork();
            proteinNetworkPre.loadRates();
            proteinNetworkPre.setLastProteinId();
        }

        if(Parameters.useIdealisedBug) {
            proteinNetworkPre.loadResponseScalingFactor();
        } if(Parameters.useRealisticBug) {
            proteinNetworkPre.loadRealisticBug();
        }

        if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            proteinNetworkPre.fillTauLogSpace();
        }

        RungeKutta rungeKuttaPre = new RungeKutta(proteinNetworkPre.getProteinIds()), rungeKuttaPost = null;
        World worldPre = new World(proteinNetworkPre.getProteinIds()), worldPost;
        World.setSensitivity();
        World.initialiseGaussians();
        if(Parameters.loadNetwork == -1) proteinNetworkPre.setSensitivity(World.getSensitivity());
        Mutagenesis mutagenesis = new Mutagenesis();
        if(Parameters.loadNetwork == -1 && !Parameters.useIdealisedBug && !Parameters.useRealisticBug) mutagenesis.fillNetwork(proteinNetworkPre);
        List<List<ArrayList<ArrayList<Float>>>> fourierCoefficients = new ArrayList<List<ArrayList<ArrayList<Float>>>>();
        List<SingleBug> singleBugs = new ArrayList<SingleBug>();

        for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
            singleBugs.add(new SingleBug());
        }

        Joiner j = Joiner.on(",");

        for(int generation = 1; generation <= Parameters.numberOfGenerations; generation++) {
            if(Parameters.printGeneration && generation % Parameters.printGenerationEvery == 0) {
                System.out.println("G" + generation);
            }

            World.updateSensitivity(generation);
            seed = getSeed();
            random = new Well44497b(seed);
            System.out.println(seed);

            try { proteinNetworkPost = (ProteinNetwork) proteinNetworkPre.clone();
            } catch(CloneNotSupportedException e) { System.exit(-47); }
            try { rungeKuttaPost = (RungeKutta) rungeKuttaPre.clone();
            } catch(CloneNotSupportedException e) { System.exit(-47); }

            mutagenesis.mutate(proteinNetworkPost, rungeKuttaPost);
            System.out.println(j.join(proteinNetworkPre.tauLogSpace));
            System.out.println(j.join(proteinNetworkPost.tauLogSpace));

            if(Parameters.useStochasticFood && Parameters.equilibrateStochasticFood) {
                World.determineEquilibrationPeriod(proteinNetworkPre.realisticBugParameters.get(3), proteinNetworkPost.realisticBugParameters.get(3));
            }

            if(Parameters.useRealisticBug) {
                worldPre.memoryCalculations(proteinNetworkPre);
                worldPre.calculateExponentials(proteinNetworkPre);
            }

            timeBefore = System.nanoTime();

            if(Parameters.curveFitting) {
                tetheredBug(worldPre, proteinNetworkPre);
            } else {
                fourierCoefficients = oneGeneration(proteinNetworkPre, rungeKuttaPre, worldPre, singleBugs, null, false);
            }

            if(Parameters.printCPUDurations) {
                System.out.println((((double) (System.nanoTime() - timeBefore)) / 1000000000));
            }

            rungeKuttaPre.resetActivityLevels();
            worldPost = new World(proteinNetworkPost.getProteinIds());

            if(Parameters.useRealisticBug) {
                worldPost.memoryCalculations(proteinNetworkPost);
                worldPost.calculateExponentials(proteinNetworkPost);
            }

            timeBefore = System.nanoTime();

            if(Parameters.curveFitting) {
                tetheredBug(worldPost, proteinNetworkPost);
            } else {
                oneGeneration(proteinNetworkPost, rungeKuttaPost, worldPost, singleBugs, fourierCoefficients, true);
            }

            if(Parameters.printCPUDurations) {
                System.out.println((((double) (System.nanoTime() - timeBefore)) / 1000000000));
                System.out.println((((double) (System.nanoTime() - timeStart)) / 1000000000));
            }

            worldPre.resetExponentials();
            worldPost.resetExponentials();
            rungeKuttaPost.resetActivityLevels();
            acceptanceProbability = World.fitnessFunction(worldPre.getFitness(), worldPost.getFitness());

            if(Parameters.printPreAndPostDurations && networkShouldBePrinted(generation)) {
                worldPre.appendToDurations(generation);
                worldPost.appendToDurations(generation);
            }

            if(acceptanceProbability == 1) {
                fitnessRejected = worldPre.getFitness();
            } else if(acceptanceProbability == 0) {
                fitnessRejected = worldPost.getFitness();
            }

            if(Parameters.printBothPreAndPost) {
                if(networkShouldBePrinted(generation)) proteinNetworkPre.appendToFile(generation);
                worldPre.appendToFile(generation, networkShouldBePrinted(generation), fitnessRejected);
                if(networkShouldBePrinted(generation)) proteinNetworkPost.appendToFile(generation);
                worldPost.appendToFile(generation, networkShouldBePrinted(generation), fitnessRejected);
            }

            if(acceptanceProbability == 1) {
                try { proteinNetworkPre = (ProteinNetwork) proteinNetworkPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }
                try { rungeKuttaPre = (RungeKutta) rungeKuttaPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }
                try { worldPre = (World) worldPost.clone();
                } catch(CloneNotSupportedException e) { System.exit(-47); }
            }

            if(!Parameters.printBothPreAndPost) {
                if(networkShouldBePrinted(generation)) proteinNetworkPre.appendToFile(generation);
                worldPre.appendToFile(generation, networkShouldBePrinted(generation), fitnessRejected);
            }

            rungeKuttaPre.resetStepSizeNextStep();

            if(Parameters.flushEvery != -1 && generation % Parameters.flushEvery == 0) {
                flush();
            }
        }

        if(Parameters.curveFitting) {
            List<List<Double>> matrix = new ArrayList<List<Double>>();
            matrix.add(worldPre.times);
            matrix.add(worldPre.effectorActivities);
            writeCSV("curve_fit.csv", matrix);
        }
    }

    private void tetheredBug(World world, ProteinNetwork proteinNetwork) {
        double stepSize = Parameters.useRealisticBug ? proteinNetwork.bugStepSize : Parameters.idealisedBugStepSize;
        Bug bug = new Bug(Parameters.driftVelocity);
        bug.initialiseFoodsPast(world.getFoodsPast(proteinNetwork, bug));
        bug.setLocalCopiesOfParameters(proteinNetwork);
        bug.equilibrateAlpha();

        while(world.getSteps() <= proteinNetwork.howManyStepsInUnitTime * Parameters.stimulusCourseDuration) {
            world.increaseSteps();
            world.setTime(world.getSteps() * stepSize);
            world.logTime();
            world.addEffectorActivity(bug.moveRealisticFasterTetheredBug(World.getInput(world.getTime() - stepSize), World.getInput(world.getTime())));
            world.advanceInput();
        }

        world.calculateFitness();
    }

    private List<List<ArrayList<ArrayList<Float>>>> oneGeneration(ProteinNetwork proteinNetwork, RungeKutta rk, World world, List<SingleBug> singleBugs, List<List<ArrayList<ArrayList<Float>>>> fourierCoefficientsToUse, boolean isPostMutagenesis) {
        RungeKutta rkClone = null;
        World worldClone = null;
        List<Thread> threads = new ArrayList<Thread>();
        List<Boolean> threadsFinished = new ArrayList<Boolean>();
        List<Double> threadsFinishedDurations = new ArrayList<Double>();
        List<List<ArrayList<ArrayList<Float>>>> fourierCoefficients = new ArrayList<List<ArrayList<ArrayList<Float>>>>();
        ArrayList<Double> finalDurations = new ArrayList<Double>();
        int finishedCount = 0;
        BitsStreamGenerator random;
        long seed;

        if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
            world.waitUntilSteadyState(rk, proteinNetwork, new Bug(-1));
        }

        if(world.isEquilibrated() || Parameters.useIdealisedBug || Parameters.useRealisticBug) {
            for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                try { rkClone = (RungeKutta) rk.clone(); } catch(CloneNotSupportedException e) {}
                try { worldClone = (World) world.clone(); } catch(CloneNotSupportedException e) {}
                seed = getSeed();
                random = new Well44497b(seed);
                worldClone.setRandom(random);

                if(Parameters.useStochasticFood && Parameters.equilibrateStochasticFood && !isPostMutagenesis) {
                    worldClone.equilibrateStochasticFood();
                }

                if(Parameters.useStochasticFood && isPostMutagenesis) {
                    worldClone.setFourierCoefficients(fourierCoefficientsToUse.get(i));
                }

                if(!Parameters.differentBugsDifferentDrifts) {
                    singleBugs.get(i).initialise(proteinNetwork, rkClone, worldClone, Parameters.driftVelocity, random);
                } else {
                    if(i < Parameters.numberOfConcurrentBugs / 2) {
                        singleBugs.get(i).initialise(proteinNetwork, rkClone, worldClone, 0, random);
                    } else {
                        singleBugs.get(i).initialise(proteinNetwork, rkClone, worldClone, Parameters.driftVelocity, random);
                    }
                }

                threads.add(new Thread(singleBugs.get(i)));
                System.out.println(threads.get(i).getId() + "," + seed);
                threads.get(i).start();
                threadsFinished.add(false);
            }

            try {
                while(finishedCount < Parameters.howManyBugsToConsider) {
                    Thread.sleep(10);
                    finishedCount = 0;

                    for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                        if(!threads.get(i).isAlive()) {
                            finishedCount++;
                            threadsFinished.set(i, true);
                        }
                    }
                }
            } catch(InterruptedException e) {}

            for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                if(threadsFinished.get(i)) {
                    threadsFinishedDurations.add(singleBugs.get(i).getWholeEstimatesCount());
                }
            }

            Collections.sort(threadsFinishedDurations);
            //System.out.println("Main thread: " + threadsFinishedDurations.get(Parameters.howManyBugsToConsider - 1) * Parameters.maxDurationFactor);

            for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                if(!threadsFinished.get(i)) {
                    singleBugs.get(i).setWholeEstimatesLimit(threadsFinishedDurations.get(Parameters.howManyBugsToConsider - 1) * Parameters.maxDurationFactor);
                }
            }

            try {
                while(finishedCount < Parameters.numberOfConcurrentBugs) {
                    Thread.sleep(10);
                    finishedCount = 0;

                    for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                        if(!threads.get(i).isAlive()) {
                            finishedCount++;

                            if(!threadsFinished.get(i)) {
                                threadsFinished.set(i, true);
                                threadsFinishedDurations.add(singleBugs.get(i).getWholeEstimatesCount());
                                Collections.sort(threadsFinishedDurations);
                                //System.out.println("Main thread: " + threadsFinishedDurations.get(Parameters.howManyBugsToConsider - 1) * Parameters.maxDurationFactor);

                                for(int j = 0; j < Parameters.numberOfConcurrentBugs; j++) {
                                    if(!threadsFinished.get(j)) {
                                        singleBugs.get(j).setWholeEstimatesLimit(threadsFinishedDurations.get(Parameters.howManyBugsToConsider - 1) * Parameters.maxDurationFactor);
                                    }
                                }
                            }
                        }
                    }
                }
            } catch(InterruptedException e) {}

            for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
                if(singleBugs.get(i).getWholeEstimatesCount() == -1) {
                    System.out.println("wholeEstimatesCount is -1, the associated thread has probably crashed, exiting simulation.");
                    System.exit(-47);
                }

                finalDurations.add(singleBugs.get(i).getWholeEstimatesCount() / Parameters.scaleFitnessBy);

                if(Parameters.useStochasticFood && !isPostMutagenesis) {
                    fourierCoefficients.add(singleBugs.get(i).getFourierCoefficients());
                }
            }

            world.setDurations(finalDurations);
            world.calculateFitness();
        } else {
            System.out.println("Steady state has not been reached - fitness will be set to 0 and the concurrent bugs will not be run.");
            world.setFitness(0);
        }

        return fourierCoefficients;
    }

    private void testSingleBugOneGeneration() {
        ProteinNetwork proteinNetwork = loadNetwork();
        proteinNetwork.loadRates();
        RungeKutta rk = new RungeKutta(proteinNetwork.getProteinIds());
        World world = new World(proteinNetwork.getProteinIds());
        World.initialiseGaussians();
        List<SingleBug> singleBugs = new ArrayList<SingleBug>();

        for(int i = 0; i < Parameters.numberOfConcurrentBugs; i++) {
            singleBugs.add(new SingleBug());
        }

        oneGeneration(proteinNetwork, rk, world, singleBugs, null, false);
        world.appendToDurations(-1);
        System.out.println("Fitness: " + world.getFitness());
    }

    private static class SingleBug implements Runnable {
        ProteinNetwork proteinNetwork;
        RungeKutta rk;
        World world;
        double wholeEstimatesCount, wholeEstimatesLimit, foodRunningSum, driftVelocity;
        BitsStreamGenerator random;

        public void initialise(ProteinNetwork proteinNetwork, RungeKutta rk, World world, double driftVelocity, BitsStreamGenerator r) {
            this.proteinNetwork = proteinNetwork;
            this.rk = rk;
            this.world = world;
            wholeEstimatesLimit = Double.MAX_VALUE;
            foodRunningSum = 0;
            this.driftVelocity = driftVelocity;
            wholeEstimatesCount = -1;
            random = r;
        }

        public void run() {
            double currentFood = 0, foodMidStep, stepSize = Parameters.useRealisticBug ? proteinNetwork.bugStepSize : Parameters.idealisedBugStepSize;
            Bug bug = new Bug(driftVelocity);
            bug.setRandom(random);

            if(Parameters.useIdealisedBug) {
                bug.initialiseFoodPast(world.stochasticFoodDistribution(bug.getPosition(), -1));
            } if(Parameters.useRealisticBug) {
                bug.initialiseFoodsPast(world.getFoodsPast(proteinNetwork, bug));
                bug.setExponentials(world.getExponentials());
                bug.setLocalCopiesOfParameters(proteinNetwork);
                bug.equilibrateAlpha();
            }

            long timeBefore = System.nanoTime();

            while(world.getTime() + rk.getStepSizeNextStep() <= wholeEstimatesLimit) {
                world.increaseSteps();

                if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
                    if(!rk.cashKarp(proteinNetwork, world, bug, false)) {
                        System.out.println("fourthOrder used!");
                        System.exit(-47);
                    }

                    world.increaseTime(rk.getStepSizeThisStep());
                } else {
                    world.setTime(world.getSteps() * stepSize);
                }

                if(world.getTime() > Parameters.maxWait) break;

                foodMidStep = world.stochasticFoodDistribution(bug.getMidPosition(), World.getIndexOfNearestStochasticTime(world.getTime() - stepSize));
                currentFood = world.stochasticFoodDistribution(bug.getPosition(), World.getIndexOfNearestStochasticTime(world.getTime()));

                if(Parameters.useIdealisedBug) {
                    bug.moveIdealised(proteinNetwork, Parameters.idealisedBugStepSize, currentFood);
                } else if(Parameters.useRealisticBug) {
                    if(!Parameters.instantaneousTumbles) {
                        bug.moveRealisticFaster(foodMidStep, currentFood);
                    } else {
                        bug.moveInstantaneousTumbles(foodMidStep, currentFood);
                    }
                } else {
                    bug.move(proteinNetwork, rk.getActivityLevel(Protein.EXTRA_INTERACTIONS + 1), rk.getStepSizeThisStep());
                }

                if(world.isTimeWhole(proteinNetwork.howManyStepsInUnitTime)) {
                    /*if(world.getTime() % 100 == 0) {
                        System.out.println(Thread.currentThread().getName() + ": " + world.getTime() + ", food: " + foodRunningSum);
                    }*/

                    if(Parameters.useStochasticFood) {
                        foodRunningSum += currentFood;
                    } else {
                        foodRunningSum += world.getFoodAt(bug.getPosition()) - Parameters.gaussianBaselineConcentration;
                    }

                    if(foodRunningSum >= Parameters.foodRunningSumLimit) {
                        wholeEstimatesCount = world.getTime();
                        //System.out.format("%s: Done! Got %f food in %f time.\n", Thread.currentThread().getName(), foodRunningSum, wholeEstimatesCount);
                        break;
                    }
                }
            }

            if(foodRunningSum < Parameters.foodRunningSumLimit) {
                wholeEstimatesCount = Integer.MAX_VALUE;
                //System.out.format("%s: Didn't finish in time (%f), only got %f.\n", Thread.currentThread().getName(), wholeEstimatesLimit, foodRunningSum);
            }

            if(Parameters.printCPUDurations) {
                if(Parameters.useIdealisedBug || Parameters.useRealisticBug) {
                    System.out.format("%d,%f,%d\n", Thread.currentThread().getId(), ((double) (System.nanoTime() - timeBefore)) / 1000000000, (int) wholeEstimatesCount);
                } else {
                    System.out.format("%d,%f,%d,%d\n", Thread.currentThread().getId(), ((double) (System.nanoTime() - timeBefore)) / 1000000000, world.getSteps(), (int) wholeEstimatesCount);
                }
            }
        }

        public double getWholeEstimatesCount() {
            return wholeEstimatesCount;
        }

        public void setWholeEstimatesLimit(double d) {
            wholeEstimatesLimit = d;
        }

        public List<ArrayList<ArrayList<Float>>> getFourierCoefficients() {
            return world.getFourierCoefficients();
        }
    }

    private ProteinNetwork loadNetwork() {
        ProteinNetwork proteinNetwork = new ProteinNetwork();
        Scanner effectors = null, receptors = null, roles = null, interactions = null, effectorMotorRates = null;

        try {
            if(Modularity.whatTodo == 3) {
                if(!Parameters.useIdealisedBug && !Parameters.useRealisticBug) {
                    effectors = new Scanner(new BufferedReader(new FileReader(filenames.get(7))));
                    receptors = new Scanner(new BufferedReader(new FileReader(filenames.get(5))));
                    roles = new Scanner(new BufferedReader(new FileReader(filenames.get(6))));
                    interactions = new Scanner(new BufferedReader(new FileReader(filenames.get(4))));
                } else if((Parameters.useIdealisedBug || Parameters.useRealisticBug) && loadGeneration != -1) {
                    effectorMotorRates = new Scanner(new BufferedReader(new FileReader(filenames.get(10)))).useDelimiter("\n");
                }
            } else if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2) {
                effectors = new Scanner(new BufferedReader(new FileReader("effectors1_" + Parameters.loadNetwork)));
                receptors = new Scanner(new BufferedReader(new FileReader("receptors1_" + Parameters.loadNetwork)));
                roles = new Scanner(new BufferedReader(new FileReader("roles1_" + Parameters.loadNetwork)));
                interactions = new Scanner(new BufferedReader(new FileReader("interactions1_" + Parameters.loadNetwork)));
            }

            String line;

            if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2 || (Modularity.whatTodo == 3 && !Parameters.useIdealisedBug && !Parameters.useRealisticBug)) {
                effectors.useDelimiter("\n");
                receptors.useDelimiter("\n");
                roles.useDelimiter("\n");
                interactions.useDelimiter("\n");
                String interactionsRegex = " ", restRegex = " = ";
                Pattern interactionsPattern = Pattern.compile(interactionsRegex), restPattern = Pattern.compile(restRegex);
                HashMap<Integer, Protein> proteins = new HashMap<Integer, Protein>();
                effectors.next();
                receptors.next();
                roles.next();
                interactions.next();
                String[] lineSplit;
                Protein protein;

                while(effectors.hasNext()) {
                    line = effectors.next();
                    if(line.equals("")) break;
                    lineSplit = restPattern.split(line);
                    protein = new Protein(Integer.parseInt(lineSplit[0]));
                    if(lineSplit[1].equals("isEffector")) protein.setEffector(true);
                    proteins.put(protein.getProteinId(), protein);
                } while(receptors.hasNext()) {
                    line = receptors.next();
                    if(line.equals("")) break;
                    lineSplit = restPattern.split(line);
                    if(lineSplit[1].equals("isReceptor")) proteins.get(Integer.parseInt(lineSplit[0])).setReceptor(true);
                } while(roles.hasNext()) {
                    line = roles.next();
                    if(line.equals("")) break;
                    lineSplit = restPattern.split(line);
                    if(lineSplit[1].equals("isPhosphatase")) proteins.get(Integer.parseInt(lineSplit[0])).switchRole();
                } while(interactions.hasNext()) {
                    line = interactions.next();
                    if(line.equals("")) break;
                    lineSplit = interactionsPattern.split(line);

                    if(lineSplit[0].equals("F")) {
                        proteins.get(Integer.parseInt(lineSplit[2])).addInteraction(Protein.RECEPTOR_SENSITIVITY, Double.parseDouble(lineSplit[4]));
                    } else if(lineSplit[0].equals(lineSplit[2])) {
                        if(lineSplit[1].equals("(deactivates)"))
                            proteins.get(Integer.parseInt(lineSplit[0])).addInteraction(0, Double.parseDouble(lineSplit[4]));
                        else
                            proteins.get(Integer.parseInt(lineSplit[0])).addInteraction(1, Double.parseDouble(lineSplit[4]));
                    } else {
                        proteins.get(Integer.parseInt(lineSplit[0])).addInteraction(Integer.parseInt(lineSplit[2]), Double.parseDouble(lineSplit[4]));
                    }
                }

                for(Protein p : proteins.values())
                    proteinNetwork.addProtein(p.getProteinId(), p);
            }

            if(Modularity.whatTodo == 3 && (Parameters.useIdealisedBug || Parameters.useRealisticBug) && loadGeneration != -1) {
                Pattern generationPattern = Pattern.compile("^" + loadGeneration + ","), commaPattern = Pattern.compile(",");
                Matcher generationMatcher;
                String[] items;

                while(effectorMotorRates.hasNext()) {
                    line = effectorMotorRates.next();
                    generationMatcher = generationPattern.matcher(line);

                    if(generationMatcher.find()) {
                        items = commaPattern.split(line);
                        Parameters.effectorMotorAssociationRate = Double.parseDouble(items[1]);
                        Parameters.effectorMotorDissociationRate = Double.parseDouble(items[2]);

                        if(Parameters.useIdealisedBug) {
                            Parameters.responseScalingFactor = Double.parseDouble(items[3]);
                        } else if(Parameters.useRealisticBug) {
                            Parameters.realisticBugConstant = Double.parseDouble(items[3]);
                            Parameters.realisticBugLinear = Double.parseDouble(items[4]);
                            Parameters.realisticBugQuadratic = Double.parseDouble(items[5]);
                            Parameters.realisticBugMemoryLength = Double.parseDouble(items[6]);
                        }

                        if(!Parameters.loadNetworkPost) {
                            break;
                        }
                    }
                }
            }
        } catch(IOException e) {
            System.out.println("Input files not present.");
            System.exit(-47);
        } finally {
            if(effectors != null) effectors.close();
            if(receptors != null) receptors.close();
            if(roles != null) roles.close();
            if(interactions != null) interactions.close();
            if(effectorMotorRates != null) effectorMotorRates.close();
        }

        return proteinNetwork;
    }

    private void loadCSV(String filename, List<List<Double>> matrix) {
        Scanner s = null;
        Pattern commaPattern = Pattern.compile(",");
        String[] items;
        String line;

        try {
            s = new Scanner(new BufferedReader(new FileReader(filename)));

            while(s.hasNext()) {
                line = s.next();
                items = commaPattern.split(line);

                if(items.length != matrix.size()) {
                    System.out.println("The supplied List<List<Double>> is of incorrect size.");
                    System.exit(-47);
                }

                for(int i = 0; i < items.length; i++) {
                    matrix.get(i).add(Double.parseDouble(items[i]));
                }
            }
        } catch(FileNotFoundException e) {
            System.out.println("File " + filename + " could not be found.");
            System.exit(-47);
        } finally {
            if(s != null) {
                s.close();
            }
        }
    }

    private void writeCSV(String filename, List<List<Double>> matrix) {
        PrintWriter csv = null;
        int noOfRowsInFirstColumn = matrix.get(0).size();
        NumberFormat f = new DecimalFormat("0.######E0");

        try {
            csv = new PrintWriter(new BufferedWriter(new FileWriter(filename)));

            for(List<Double> l : matrix) {
                if(l.size() != noOfRowsInFirstColumn) {
                    System.out.println("Unequal numbers of rows in columns of the matrix.");
                    System.exit(-47);
                }
            }

            for(int row = 0; row < noOfRowsInFirstColumn; row++) {
                for(int column = 0; column < matrix.size(); column++) {
                    if(column == matrix.size() - 1) {
                        csv.println(f.format(matrix.get(column).get(row)));
                    } else {
                        csv.print(f.format(matrix.get(column).get(row)) + ",");
                    }
                }
            }
        } catch(IOException e) {
            System.out.println("File " + filename + " could not be written.");
            System.exit(-47);
        } finally {
            if(csv != null) {
                csv.close();
            }
        }
    }

    private void initialiseFilenames(int jobNumber, int parameterSet, int run) {
        filenames.add(0, "output");
        filenames.add(1, "fitnesses");
        filenames.add(2, "orkun_fitness_components");
        filenames.add(3, "network");
        filenames.add(4, "interactions");
        filenames.add(5, "receptors");
        filenames.add(6, "roles");
        filenames.add(7, "effectors");
        filenames.add(8, "fitness_rejected");
        filenames.add(9, "fixation_probabilities");
        filenames.add(10, "effector_motor_rates");
        filenames.add(11, "durations");
        filenames.add(12, "positions");
        filenames.add(13, "tumbling_probabilities");
        filenames.add(14, "run_durations");
        filenames.add(15, "tumble_durations");
        filenames.add(16, "alphas");
        filenames.add(17, "running_probabilities");
        filenames.add(18, "convolutions");
        int filenamesSize = filenames.size();

        if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2) {
            for(int i = 0; i < filenamesSize; i++) {
                filenames.set(i, filenames.get(i) + parameterSet + "_" + run);
            }
        } else if(Modularity.whatTodo == 3 || Modularity.whatTodo == 9) {
            for(int i = 0; i < filenamesSize; i++) {
                filenames.set(i, filenames.get(i) + parameterSet + "_" + jobNumber);
            }
        }

        if(Modularity.whatTodo == 1) {
            for(int i = 0; i < filenamesSize; i++) {
                filenames.set(i, "./" + jobNumber + "/" + filenames.get(i));
            }
        }
    }

    private void clearFiles() {
        File file;

        for(int i = 0; i < filenames.size(); i++) {
            if(((i == 1 || (i >= 3 && i <= 9)) && (Modularity.whatTodo == 9 || Modularity.whatTodo == 3))
                    || (i == 0 && (Modularity.whatTodo == 9 || (Parameters.whatToSimulate == 2 && (Modularity.whatTodo == 1 || Modularity.whatTodo == 2))))
                    || (i == 2 && (Modularity.whatTodo == 9 || Parameters.whatToSimulate == 2 || !Parameters.fitnessCalculation[3]))
                    || (i == 10 && (Modularity.whatTodo == 9 || Parameters.whatToSimulate == 1 || (Modularity.whatTodo == 3 && Parameters.whatToSimulate == 2)))
                    || (i == 11 && (Parameters.whatToSimulate == 1 || (Modularity.whatTodo == 3 && Parameters.whatToSimulate == 2)))
                    || ((i >= 12 && i <= 18) && (Modularity.whatTodo == 9 || Parameters.whatToSimulate == 1 || ((Modularity.whatTodo == 1 || Modularity.whatTodo == 2) && Parameters.whatToSimulate == 2)))) {
                continue;
            }

            file = new File(filenames.get(i));
            file.delete();
            try { file.createNewFile(); } catch(IOException e) {}
        }
    }

    private void initialiseFiles() {
        try {
            if(Modularity.whatTodo != 9 && !(Parameters.whatToSimulate == 2 && (Modularity.whatTodo == 1 || Modularity.whatTodo == 2))) {
                output = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(0), true)));
            } if((Modularity.whatTodo == 1 || Modularity.whatTodo == 2) && Parameters.whatToSimulate == 1) {
                output.format("\"Generation\",\"Time\",\"Output\"\n");
            } if(Modularity.whatTodo == 3) {
                output.format("\"Time\",\"Input\"");
            }

            if(Parameters.whatToSimulate == 1 && Parameters.fitnessCalculation[3] && Modularity.whatTodo != 9) {
                orkunComponents = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(2), true)));

                if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2) {
                    orkunComponents.format("\"Generation\",\"Time\",\"Inverted\",\"Adaptive\",\"Constant\"\n");
                } else if(Modularity.whatTodo == 3) {
                    orkunComponents.format("\"Time\",\"Inverted\",\"Adaptive\",\"Constant\"\n");
                }
            }

            if(Parameters.whatToSimulate == 2 && Modularity.whatTodo != 3) {
                durations = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(11), true)));
            }

            if(Parameters.whatToSimulate == 2 && Modularity.whatTodo == 3) {
                positions = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(12), true)));
                positions.format("\"Time\",\"Position\"\n");
                tumblingProbabilities = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(13))));
                tumblingProbabilities.format("\"Time\",\"Probability\"\n");
                runDurations = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(14))));
                runDurations.format("\"Duration\"\n");
                tumbleDurations = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(15))));
                tumbleDurations.format("\"Duration\"\n");
                alphas = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(16))));
                alphas.format("\"Time\",\"Alpha\"\n");
                runningProbabilities = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(17))));
                runningProbabilities.format("\"Time\",\"Probability\"\n");
                convolutions = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(18))));
                convolutions.format("\"Time\",\"Integral\"\n");
            }

            if(Modularity.whatTodo == 1 || Modularity.whatTodo == 2) {
                fitnesses = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(1), true)));

                if(Parameters.fitnessCalculation[2])
                    fitnesses.format("\"Generation\",\"Linear\",\"Adaptive\",\"Flat\"\n");
                else
                    fitnesses.format("\"Generation\",\"Fitness\"\n");

                network = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(3), true)));
                interactions = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(4), true)));
                receptors = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(5), true)));
                roles = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(6), true)));
                effectors = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(7), true)));
                fitnessRejectedWriter = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(8), true)));
                fitnessRejectedWriter.format("\"Generation\",\"Fitness\"\n");
                fixationProbabilities = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(9), true)));
                fixationProbabilities.format("\"Generation\",\"Probability\"\n");

                if(Parameters.whatToSimulate == 2) {
                    effectorMotorRates = new PrintWriter(new BufferedWriter(new FileWriter(filenames.get(10), true)));

                    if(Parameters.useIdealisedBug) {
                        effectorMotorRates.format("\"Generation\",\"Association_Rate\",\"Dissociation_Rate\",\"Response_Scaling_Factor\"\n");
                    } else if(Parameters.useRealisticBug) {
                        effectorMotorRates.format("\"Generation\",\"Association_Rate\",\"Dissociation_Rate\",\"Constant\",\"Linear\",\"Quadratic\",\"Memory_Length\",\"epsilon\"\n");
                    } else {
                        effectorMotorRates.format("\"Generation\",\"Association_Rate\",\"Dissociation_Rate\"\n");
                    }
                }
            }

            urandom = new DataInputStream(new FileInputStream("/dev/urandom"));
        } catch(IOException e) {}
    }

    private void closeFiles() {
        if(fitnesses != null) {
            fitnesses.close();
        } if(output != null) {
            output.close();
        } if(orkunComponents != null) {
            orkunComponents.close();
        } if(network != null) {
            network.close();
        } if(interactions != null) {
            interactions.close();
        } if(receptors != null) {
            receptors.close();
        } if(roles != null) {
            roles.close();
        } if(effectors != null) {
            effectors.close();
        } if(fitnessRejectedWriter != null) {
            fitnessRejectedWriter.close();
        } if(fixationProbabilities != null) {
            fixationProbabilities.close();
        } if(effectorMotorRates != null) {
            effectorMotorRates.close();
        } if(durations != null) {
            durations.close();
        } if(positions != null) {
            positions.close();
        } if(tumblingProbabilities != null) {
            tumblingProbabilities.close();
        } if(runDurations != null) {
            runDurations.close();
        } if(tumbleDurations != null) {
            tumbleDurations.close();
        } if(alphas != null) {
            alphas.close();
        } if(runningProbabilities != null) {
            runningProbabilities.close();
        } if(convolutions != null) {
            convolutions.close();
        } if(urandom != null) {
            try {
                urandom.close();
            } catch(IOException e) {}
        }
    }

    private void testRKExponentialDecay() {
        Plotting plot = new Plotting();
        plot.initialise();
        RungeKutta rk = new RungeKutta(null);
        ProteinNetwork pn = new ProteinNetwork();
        World world = new World(pn.getProteinIds());

        System.out.println("beginning\n");

        while(world.getTime() <= 150) {
            if(!rk.cashKarp(pn, world.getInput())) {
                rk.fourthOrder(pn, world.getInput());
                System.out.print("fourthOrder(): ");
            }
            else {
                System.out.print("cashKarp(): ");
            }

            world.increaseTime(rk.getStepSizeThisStep());
            System.out.println("time " + world.getTime() + ", h " + rk.getStepSizeThisStep() + ", goal " + Math.exp(-0.2 * world.getTime()) + ", output " + rk.getActivityLevel(0));
            plot.updateTrace(0, world.getTime(), rk.getActivityLevel(0));
            plot.updateTrace(3, world.getTime(), Math.exp(-0.2 * world.getTime()));
        }

        System.out.println("end\n");
    }

    private void testStochasticFood() {
        Set<Integer> set = new HashSet<Integer>();
        set.add(4);
        World world = new World(set);
        world.setRandom(Modularity.random);
        List<Double> foodDistribution = new ArrayList<Double>();
        double stepSize = 1 / Parameters.foodEstimatesPerUnitWorldLength, food, foodSum = 0;

        for(double d = 0; d <= Parameters.worldLength; d += stepSize) {
            foodDistribution.add(world.stochasticFoodDistribution(d, -1));
        }

        PlottingStochastic plot = new PlottingStochastic(0, foodDistribution);
        JFrame jFrame = new JFrame("Chemotaxis");
        jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        jFrame.add(plot);
        jFrame.pack();
        jFrame.setVisible(true);

        for(int i = 0; i < 100000; i++) {
            if(i % 1000 == 0) System.out.println(i);
            world.nextFourierCoefficient();
            foodDistribution.clear();

            for(double d = 0; d <= Parameters.worldLength; d += stepSize) {
                food = world.stochasticFoodDistribution(d, -1);
                foodDistribution.add(food);
                foodSum += food;
            }

            plot.moveBug(0);
            plot.moveBugRepaint();

            /*try {
                Thread.sleep(10);
            } catch(InterruptedException e) {}*/
        }

        //why 1000000? because time = 100000 and spatial spacing = 10. in other words, this is the average integral of food in the space.
        System.out.println(foodSum / 1000000);
        jFrame.dispose();
    }

    private void testRKProteinNetwork() {
        ProteinNetwork pn = new ProteinNetwork();
        pn.initialise();
        /*pn.getProtein(3).switchRole();
        pn.getProtein(3).addInteraction(0, 1);
        pn.getProtein(3).addInteraction(5, 1);*/
        pn.getProtein(3).addInteraction(2, 1);
        pn.getProtein(3).addInteraction(4, 1);
        pn.getProtein(4).addInteraction(5, 1);
        pn.print();

        RungeKutta rk = new RungeKutta(null);
        /*rk.setActivityLevel(3, 1);
        rk.setActivityLevel(5, 1);*/
        rk.setActivityLevel(3, 0);
        rk.setActivityLevel(5, 0);

        Dynamics d = new Dynamics(pn, 1);
        double[] activityLevels = new double[4];
        /*activityLevels[0] = 1;
        activityLevels[1] = 0;
        activityLevels[2] = 1;
        activityLevels[3] = 0;*/
        activityLevels[0] = 0;
        activityLevels[1] = 0;
        activityLevels[2] = 0;
        activityLevels[3] = 0;
        FlanagansRungeKutta.H = 0.0001;

        for(int i = 0; i < 50; i++) {
            rk.cashKarp(pn, 1);
            //rk.fourthOrder(pn, 1);
            System.out.println(rk.getActivityLevel(3) + ", " + rk.getActivityLevel(4) + ", " + rk.getActivityLevel(5) + ", " + rk.getActivityLevel(6));

            activityLevels = FlanagansRungeKutta.cashKarp(d, 0, activityLevels, 0.0001, FlanagansRungeKutta.H, Parameters.absoluteTolerance, Parameters.relativeTolerance, Parameters.maxNumberOfIterations);
            System.out.println("FlanagansRungeKutta\nH " + FlanagansRungeKutta.H);
            //activityLevels = FlanagansRungeKutta.fourthOrder(d, 0, activityLevels, 0.0001, 0.0001);
            System.out.println(activityLevels[0] + ", " + activityLevels[1] + ", " + activityLevels[2] + ", " + activityLevels[3]);

            System.out.println();
        }
    }

    private void testMutagenesis() {
        Mutagenesis m = new Mutagenesis();
        ProteinNetwork pn = new ProteinNetwork();
        pn.initialise();
        RungeKutta rk = new RungeKutta(null);

        System.out.println("beginning\n");

        for(int i = 0; i < 25000; i++) {
            m.mutate(pn, rk);
            pn.print();
        }

        System.out.println("end\n");
    }

    private void test() {
        ProteinNetwork pn = new ProteinNetwork();
        pn.initialise();
        pn.setSensitivity(1);
        /*pn.getProtein(3).addInteraction(1, 3.76);
        pn.getProtein(3).addInteraction(4, 0.08);
        pn.getProtein(4).addInteraction(3, 0.17);*/
        pn.getProtein(3).addInteraction(4, 1);
        pn.print();
        RungeKutta rk = new RungeKutta(null);

        for(int i = 0; i < 50; i++) {
            rk.cashKarp(pn, 0);
            System.out.println(i + ", " + rk.getActivityLevel(3) + ", " + rk.getActivityLevel(4) + ", " + rk.getActivityLevel(5) + ", " + rk.getActivityLevel(6));
        }
        for(int i = 0; i < 50; i++) {
            rk.cashKarp(pn, 1);
            System.out.println(i + ", " + rk.getActivityLevel(3) + ", " + rk.getActivityLevel(4) + ", " + rk.getActivityLevel(5) + ", " + rk.getActivityLevel(6));
        }
        for(int i = 0; i < 50; i++) {
            rk.cashKarp(pn, 0);
            System.out.println(i + ", " + rk.getActivityLevel(3) + ", " + rk.getActivityLevel(4) + ", " + rk.getActivityLevel(5) + ", " + rk.getActivityLevel(6));
        }
    }

    private long getSeed() {
        long seed = 0;

        try {
            seed = urandom.readLong();
        } catch(IOException e) {
            System.out.println("There was a problem getting a seed from /dev/urandom.");
        }

        return seed;
    }
}
