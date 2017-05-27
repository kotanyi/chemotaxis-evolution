import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.BitsStreamGenerator;
import java.util.*;

public class Bug {
    private double position, velocity, driftVelocity, foodPast, alphaNought, beta, alpha, odd, oddA, even, evenA, twiceStepSize, a, aSimpson, exp2, exp2Modified, exp1ModifiedA, exp1ModifiedAB, simpsonOdd, simpsonEven, midPosition, positionPast, convolution;
    private boolean isTumbling;
    public static final double VELOCITY_MAX = 1;
    private Queue<Double> foodsPast;
    private List<Double> exponentials;
    private BitsStreamGenerator random;

    public Bug(double driftVelocity) {
        position = 0;
        velocity = VELOCITY_MAX;
        this.driftVelocity = driftVelocity;
        isTumbling = false;
        odd = oddA = even = evenA = 0;
        midPosition = position;
    }

    public void initialiseFoodPast(double foodPast) {
        if(Parameters.equilibrateStochasticFood) {
            this.foodPast = foodPast;
        } else {
            this.foodPast = Parameters.stochasticBaseline;
        }
    }

    public void initialiseFoodsPast(Queue<Double> foodsPast) {
        this.foodsPast = foodsPast;
    }

    private void checkBoundaries() {
        if(position >= Parameters.worldLength) {
            position -= Parameters.worldLength;
        } else if(position < 0) {
            position += Parameters.worldLength;
        }

        if(midPosition >= Parameters.worldLength) {
            midPosition -= Parameters.worldLength;
        } else if(midPosition < 0) {
            midPosition += Parameters.worldLength;
        }
    }

    public double getPosition() {
        return position;
    }

    public boolean isTumbling() {
        return isTumbling;
    }

    public double getVelocity() {
        return velocity;
    }

    public double getDriftVelocity() {
        return driftVelocity;
    }

    public void move(ProteinNetwork pn, double effectorProteinActivity, double stepSize) {
        double alphaEffector = pn.effectorMotorRates.get(0) * effectorProteinActivity;
        double alphaEffectorPlusBeta = alphaEffector + pn.effectorMotorRates.get(1);
        double probabilityTumbling;

        if(alphaEffectorPlusBeta == 0) {
            probabilityTumbling = 0;
        } else {
            probabilityTumbling = (alphaEffector / alphaEffectorPlusBeta) * (1 - Math.exp(-alphaEffectorPlusBeta * stepSize));
        }

        double random = this.random.nextDouble();

        if(random < probabilityTumbling) {
            isTumbling = true;
        } else {
            if(isTumbling) {
                isTumbling = false;
                setNewVelocity();
            }
        }

        if(isTumbling) {
            position += driftVelocity * stepSize;
        } else {
            position += (velocity + driftVelocity) * stepSize;
        }

        checkBoundaries();
    }

    public double moveIdealised(ProteinNetwork pn, double stepSize, double food) {
        double alpha = -1;

        if(Parameters.responseFunction[0]) {
            alpha = pn.effectorMotorRates.get(0) + pn.responseScalingFactor * food;
        } else if(Parameters.responseFunction[1]) {
            alpha = pn.effectorMotorRates.get(0) - pn.responseScalingFactor * (food - foodPast) / stepSize;
            foodPast = food;
        }

        double alphaPlusBeta = alpha + pn.effectorMotorRates.get(1);
        double probabilityTumbling;

        if(alphaPlusBeta == 0) {
            probabilityTumbling = 0;
        } else {
            probabilityTumbling = (alpha / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * stepSize));
        }

        double random = this.random.nextDouble();

        if(random < probabilityTumbling) {
            isTumbling = true;
        } else {
            if(isTumbling) {
                isTumbling = false;
                setNewVelocity();
            }
        }

        if(isTumbling) {
            position += driftVelocity * stepSize;
        } else {
            position += (velocity + driftVelocity) * stepSize;
        }

        checkBoundaries();
        return probabilityTumbling;
    }

    public void setExponentials(List<Double> l) {
        exponentials = l;
    }

    public void setLocalCopiesOfParameters(ProteinNetwork pn) {
        alphaNought = pn.effectorMotorRates.get(0);
        beta = pn.effectorMotorRates.get(1);
        twiceStepSize = 2 * pn.bugStepSize;
        double tau = pn.realisticBugParameters.get(3), b = pn.realisticBugParameters.get(1) / (tau * tau), exp1 = Math.exp(-pn.bugStepSize / tau);
        a = pn.realisticBugParameters.get(0) / tau;
        aSimpson = a / 3.0 * pn.bugStepSize;
        exp1ModifiedA = a * exp1;
        exp1ModifiedAB = (a + b * pn.bugStepSize) * exp1;
        exp2 = Math.exp(-2.0 * pn.bugStepSize / tau);
        exp2Modified = 2.0 * b * pn.bugStepSize / a * exp2;
        simpsonOdd = 4.0 / 3.0 * pn.bugStepSize;
        simpsonEven = 2.0 / 3.0 * pn.bugStepSize;
    }

    public double moveRealistic(ProteinNetwork pn, double f1, double f2) {
        alpha = 0;
        foodsPast.offer(f1);
        foodsPast.offer(f2);
        foodsPast.poll();
        foodsPast.poll();
        int i = 0;

        for(Double d : foodsPast) {
            alpha += exponentials.get(i) * d;
            i++;
        }

        alpha = alpha * pn.bugStepSize / 3 + alphaNought;

        if(alpha < 0)
            alpha = 0;

        double alphaPlusBeta = alpha + beta;
        double probabilityTumbling;

        if(alphaPlusBeta == 0) {
            probabilityTumbling = 0;
        } else {
            if(!isTumbling)
                probabilityTumbling = (alpha / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * twiceStepSize));
            else
                probabilityTumbling = (beta / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * twiceStepSize));
        }

        double random = this.random.nextDouble();

        if(isTumbling && random < probabilityTumbling) {
            isTumbling = false;
            setNewVelocity();
        } else if(!isTumbling && random < probabilityTumbling) {
            isTumbling = true;
        }

        positionPast = position;

        if(isTumbling) {
            position += driftVelocity * twiceStepSize;
        } else {
            position += (velocity + driftVelocity) * twiceStepSize;
        }

        midPosition = (positionPast + position) / 2.0;
        checkBoundaries();
        return probabilityTumbling;
    }

    public void calculateOddEven(double f1, double f2) {
        odd = exp2 * odd + exp2Modified * oddA + exp1ModifiedAB * f1;
        even = exp2 * even + exp2Modified * evenA + a * f2;
        oddA = exp2 * oddA + exp1ModifiedA * f1;
        evenA = exp2 * evenA + a * f2;
    }

    public void equilibrateAlpha() {
        List<Double> foodsPastList = new ArrayList<Double>(foodsPast);
        int max = foodsPastList.size();

        for(int i = 2; i < max; i += 2) {
            calculateOddEven(foodsPastList.get(i - 1), foodsPastList.get(i));
        }
    }

    public double moveRealisticFaster(double f1, double f2) {
        calculateOddEven(f1, f2);
        convolution = simpsonOdd * odd + simpsonEven * even - aSimpson * f2;
        alpha = convolution + alphaNought;

        if(alpha < 0) {
            alpha = 0;
        }

        double probability, alphaPlusBeta = alpha + beta, random = this.random.nextDouble();

        if(isTumbling) {
            probability = (beta / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * twiceStepSize));

            if(random < probability) {
                isTumbling = false;
                setNewVelocity();
                moveRealisticFasterMove();
            }
        } else {
            probability = (alpha / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * twiceStepSize));

            if(random < probability) {
                isTumbling = true;
                midPosition = position;
            } else {
                moveRealisticFasterMove();
            }
        }

        return probability;
    }

    public double moveRealisticFasterTetheredBug(double f1, double f2) {
        calculateOddEven(f1, f2);
        alpha = simpsonOdd * odd + simpsonEven * even - aSimpson * f2 + alphaNought;

        if(alpha < 0) {
            alpha = 0;
        }

        double alphaPlusBeta = alpha + beta;
        return (alpha / alphaPlusBeta) * (1 - Math.exp(-alphaPlusBeta * twiceStepSize));
    }

    private void moveRealisticFasterMove() {
        positionPast = position;
        position += velocity * twiceStepSize;
        midPosition = (positionPast + position) / 2.0;
        checkBoundaries();
    }

    public double moveInstantaneousTumbles(double f1, double f2) {
        calculateOddEven(f1, f2);
        convolution = simpsonOdd * odd + simpsonEven * even - aSimpson * f2;
        alpha = convolution + alphaNought;

        if(alpha <= 0) {
            moveRealisticFasterMove();
            return 0;
        } else {
            AbstractRealDistribution exp = new ExponentialDistribution(random, 1/alpha, ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
            List<Double> runDurations = new ArrayList<Double>();
            double sumDurations = 0, runDuration, oddSum = 0, evenSum = 0, stepSize = twiceStepSize / 2;

            while(true) {
                runDuration = exp.sample();
                sumDurations += runDuration;

                if(sumDurations < twiceStepSize) {
                    runDurations.add(runDuration);
                } else {
                    runDurations.add(twiceStepSize - (sumDurations - runDuration));
                    break;
                }
            }

            if(runDurations.size() == 1) {
                moveRealisticFasterMove();
            } else {
                //position
                for(int i = 0; i < runDurations.size(); i++) {
                    if(i % 2 == 0) {
                        evenSum += runDurations.get(i);
                    } else {
                        oddSum += runDurations.get(i);
                    }
                }

                positionPast = position;
                position += velocity * (evenSum - oddSum);

                //midPosition
                sumDurations = evenSum = oddSum = 0;

                for(int i = 0; i < runDurations.size(); i++) {
                    runDuration = runDurations.get(i);

                    if(i % 2 == 0) {
                        evenSum += runDuration;
                    } else {
                        oddSum += runDuration;
                    }

                    if(stepSize >= sumDurations && stepSize < sumDurations + runDuration) {
                        if(i % 2 == 0) {
                            evenSum -= sumDurations + runDuration - stepSize;
                        } else {
                            oddSum -= sumDurations + runDuration - stepSize;
                        }

                        break;
                    }

                    sumDurations += runDuration;
                }

                midPosition = positionPast + velocity * (evenSum - oddSum);
                checkBoundaries();

                if(runDurations.size() % 2 == 0) {
                    velocity *= -1;
                }
            }

            return alpha;
        }
    }

    public double getMidPosition() {
        return midPosition;
    }

    public void nowWhatRandom() {
        double random = Modularity.random.nextDouble();

        if(!isTumbling) {
            if(random < 0.2) {
                isTumbling = true;
            } else {
                move();
            }
        } else {
            if(random < 0.1) {
                isTumbling = false;
                setNewVelocity();
                move();
            }
        }
    }

    private void move() {
        position += velocity;
        checkBoundaries();
    }

    private void setNewVelocity() {
        if(Parameters.pretend2D) {
            velocity = Math.cos(2 * Math.PI * this.random.nextDouble()) * VELOCITY_MAX;
        } else {
            velocity = this.random.nextInt(2) * 2 - 1;
        }
    }

    public void setRandom(BitsStreamGenerator r) {
        random = r;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getConvolution() {
        return convolution;
    }
}

