import java.util.HashSet;
import java.util.Set;

public class Mutagenesis {
    public void mutate(ProteinNetwork proteinNetwork, RungeKutta rungeKutta) {
        boolean mutagenesisSuccessful = false;
        Parameters.normaliseMutationRates(proteinNetwork.size(), proteinNetwork.getReceptorCount());

        while(!mutagenesisSuccessful) {
            double mutation = Modularity.random.nextDouble();

            if(mutation >= 0 && mutation < Parameters.thresholds.get(0)) {
                if(mutateRandomInteraction(proteinNetwork)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(0) && mutation < Parameters.thresholds.get(1)) {
                if(mutateRandomSelfInteraction(proteinNetwork)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(1) && mutation < Parameters.thresholds.get(2)) {
                if(mutateSensitivity(proteinNetwork)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(2) && mutation < Parameters.thresholds.get(3)) {
                if(duplicateRandomProtein(proteinNetwork, rungeKutta)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(3) && mutation < Parameters.thresholds.get(4)) {
                if(deleteRandomProtein(proteinNetwork, rungeKutta)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(4) && mutation < Parameters.thresholds.get(5)) {
                switchRole(proteinNetwork);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(5) && mutation < Parameters.thresholds.get(6)) {
                if(newInteractionByRecruitment(proteinNetwork, rungeKutta)) mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(6) && mutation < Parameters.thresholds.get(7)) {
                mutateAssociationRate(proteinNetwork);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(7) && mutation < Parameters.thresholds.get(8)) {
                mutateDissociationRate(proteinNetwork);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(8) && mutation < Parameters.thresholds.get(9)) {
                mutateResponseScalingFactor(proteinNetwork);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(9) && mutation < Parameters.thresholds.get(10)) {
                mutateRealisticBugParameter(proteinNetwork, 0);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(10) && mutation < Parameters.thresholds.get(11)) {
                mutateRealisticBugParameter(proteinNetwork, 1);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(11) && mutation < Parameters.thresholds.get(12)) {
                mutateRealisticBugParameter(proteinNetwork, 2);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(12) && mutation < Parameters.thresholds.get(13)) {
                mutateRealisticBugParameter(proteinNetwork, 3);
                mutagenesisSuccessful = true;
            } else if(mutation >= Parameters.thresholds.get(13) && mutation <= 1) {
                mutateEpsilon(proteinNetwork);
                mutagenesisSuccessful = true;
            }
        }
    }

    private boolean mutateRandomInteraction(ProteinNetwork proteinNetwork) {
        int from = proteinNetwork.getRandomProteinKey(null, null);
        Set<Integer> exclude = new HashSet<Integer>();
        exclude.add(from);
        int to = proteinNetwork.getRandomProteinKey(null, exclude);
        if(to == -2) return false;
        return mutateRandomInteraction(proteinNetwork, from, to);
    }

    private boolean mutateRandomSelfInteraction(ProteinNetwork proteinNetwork) {
        return mutateRandomInteraction(proteinNetwork, proteinNetwork.getRandomProteinKey(null, null), Modularity.random.nextInt(2));
    }

    private boolean mutateSensitivity(ProteinNetwork proteinNetwork) {
        int to = proteinNetwork.getRandomProteinKey(proteinNetwork.getReceptorIds(), null);
        return mutateRandomInteraction(proteinNetwork, to, Protein.RECEPTOR_SENSITIVITY);
    }

    private boolean duplicateRandomProtein(ProteinNetwork proteinNetwork, RungeKutta rungeKutta) {
        if(proteinNetwork.size() == Parameters.maxNumberOfProteins)
            return false;

        HashSet<Integer> exclude = new HashSet<Integer>();
        exclude.add(Protein.EXTRA_INTERACTIONS + 1); //Protein.EXTRA_INTERACTIONS + 1 = effector; if more effectors, can use exclude.addAll(proteinNetwork.getEffectors());
        int proteinId = proteinNetwork.getRandomProteinKey(null, exclude);

        if(proteinId == -2) { //code for "after i excluded all the proteins, there were no other proteins left and i therefore could not pick one"
            return false;
        }
        else {
            Protein duplicate = null;
            try { duplicate = (Protein) proteinNetwork.getProtein(proteinId).clone(); }
            catch(CloneNotSupportedException e) { System.exit(-47); }
            int duplicateId = proteinNetwork.increaseLastProteinId();
            duplicate.setProteinId(duplicateId);
            proteinNetwork.addProtein(duplicateId, duplicate);
            proteinNetwork.updateInteractionsAfterDuplication(proteinId, duplicateId);
            rungeKutta.addProtein(duplicateId);
            //System.out.println("Protein duplicated: " + proteinId + ", duplicateId: " + duplicateId);
            return true;
        }
    }

    private boolean deleteRandomProtein(ProteinNetwork proteinNetwork, RungeKutta rungeKutta) {
        if(proteinNetwork.size() <= Parameters.minNumberOfProteins)
            return false;

        Set<Integer> exclude = new HashSet<Integer>();
        exclude.add(Protein.EXTRA_INTERACTIONS + 1); //Protein.EXTRA_INTERACTIONS + 1 = effector; if more effectors, can use exclude.addAll(proteinNetwork.getEffectors());
        if(proteinNetwork.getReceptorCount() == 1) exclude.add(proteinNetwork.getReceptorId());
        int proteinId = proteinNetwork.getRandomProteinKey(null, exclude);

        if(proteinId == -2) { //code for "after i excluded all the proteins, there were no other proteins left and i therefore could not pick one"
            return false;
        }
        else {
            proteinNetwork.removeProtein(proteinId);

            if(proteinNetwork.size() < Parameters.minNumberOfProteins) {
                System.out.println("Protein network size: " + proteinNetwork.size());
            }

            proteinNetwork.updateInteractionsAfterDeletion(proteinId);
            rungeKutta.deleteProtein(proteinId);
            //System.out.println("Protein deleted: " + proteinId);
            return true;
        }
    }

    private void switchRole(ProteinNetwork proteinNetwork) {
        int proteinId = proteinNetwork.getRandomProteinKey(null, null);
        proteinNetwork.getProtein(proteinId).switchRole();
        //System.out.println("Protein switched: " + proteinId);
    }

    private boolean newInteractionByRecruitment(ProteinNetwork proteinNetwork, RungeKutta rungeKutta) {
        if(proteinNetwork.size() == Parameters.maxNumberOfProteins)
            return false;

        int proteinId = proteinNetwork.increaseLastProteinId();
        boolean isKinase = Modularity.random.nextBoolean();
        Protein protein = new Protein(proteinId, isKinase);
        proteinNetwork.addProtein(proteinId, protein);
        rungeKutta.addProtein(proteinId);
        boolean initiatesInteraction = Modularity.random.nextBoolean();
        Set<Integer> exclude = new HashSet<Integer>();
        exclude.add(proteinId);
        int otherProtein = proteinNetwork.getRandomProteinKey(null, exclude);
        if(otherProtein == -2) return false;
        double mutation;
        do { mutation = Modularity.random.nextDouble(); } while(mutation == 0);

        if(initiatesInteraction) {
            //System.out.println("from " + proteinId + ", to " + otherProtein + ", mutation " + mutation);
            proteinNetwork.getProtein(proteinId).addInteraction(otherProtein, mutation);
        } else {
            //System.out.println("from " + otherProtein + ", to " + proteinId + ", mutation " + mutation);
            proteinNetwork.getProtein(otherProtein).addInteraction(proteinId, mutation);
        }

        return true;
    }

    private boolean mutateRandomInteraction(ProteinNetwork proteinNetwork, int from, int to) {
        double mutation;

        if(proteinNetwork.getProtein(from).containsInteraction(to)) {
            do {
                mutation = proteinNetwork.getProtein(from).getInteraction(to) * (Parameters.deletionToAdditionRatio * Modularity.random.nextDouble() - (Parameters.deletionToAdditionRatio - 1));
            } while(mutation == 0);

            //System.out.println("from " + from + ", to " + to + ", mutation " + mutation);

            if(proteinNetwork.getProtein(from).getInteraction(to) + mutation <= 0) {
                if(Parameters.allowChangeOfTopology) {
                    proteinNetwork.getProtein(from).deleteInteraction(to);
                    return true; //the return true statements could be removed and one return true statement could be put at the end of the method. however, if deleteInteraction took place, the top-most else would also execute.
                } else {
                    return false;
                }
            }
            if(proteinNetwork.getProtein(from).getInteraction(to) + mutation >= Parameters.maxInteractionStrength) {
                if(proteinNetwork.getProtein(from).getInteraction(to) == Parameters.maxInteractionStrength) {
                    return false;
                } else {
                    proteinNetwork.getProtein(from).addInteraction(to, Parameters.maxInteractionStrength);
                    return true;
                }
            }
            else {
                proteinNetwork.getProtein(from).addInteraction(to, proteinNetwork.getProtein(from).getInteraction(to) + mutation);
                return true;
            }
        }
        else {
            if(Parameters.allowChangeOfTopology) {
                do {
                    mutation = Parameters.deletionToAdditionRatio * Modularity.random.nextDouble() - (Parameters.deletionToAdditionRatio - 1);
                } while(mutation == 0);

                if(mutation < 0) {
                    return false;
                } else {
                    //System.out.println("from " + from + ", to " + to + ", mutation " + mutation);
                    proteinNetwork.getProtein(from).addInteraction(to, mutation);
                    return true;
                }
            } else {
                return false;
            }
        }
    }

    private void mutateAssociationRate(ProteinNetwork pn) {
        mutateLogSpace(pn, 0, 0);
        pn.effectorMotorRates.set(0, Math.exp(pn.tauLogSpace.get(0)));
    }

    private void mutateDissociationRate(ProteinNetwork pn) {
        mutateLogSpace(pn, 1, 0);
        pn.effectorMotorRates.set(1, Math.exp(pn.tauLogSpace.get(1)));
    }

    private void mutateResponseScalingFactor(ProteinNetwork pn) {
        mutateLogSpace(pn, 2, 0);
        pn.responseScalingFactor = Math.exp(pn.tauLogSpace.get(2));
    }

    private void mutateRealisticBugParameter(ProteinNetwork pn, int index) {
        if(index == 0) {
            if(!Parameters.useNewWayOfMutatingAAndB) {
                mutateLogSpace(pn, 3, 0);
                pn.realisticBugParameters.set(index, Math.exp(pn.tauLogSpace.get(3)));
            } else {
                mutateLogSpaceAllowNegative(pn, 3);
                pn.realisticBugParameters.set(index, Math.signum(pn.tauLogSpace.get(3)) * Math.exp(Math.abs(pn.tauLogSpace.get(3)) + Parameters.minimumPossibleDecimalValue));
            }
        } else if(index == 1) {
            if(!Parameters.useNewWayOfMutatingAAndB) {
                mutateLogSpace(pn, 4, 0);
                pn.realisticBugParameters.set(index, Math.exp(pn.tauLogSpace.get(4)));
            } else {
                mutateLogSpaceAllowNegative(pn, 4);
                pn.realisticBugParameters.set(index, Math.signum(pn.tauLogSpace.get(4)) * Math.exp(Math.abs(pn.tauLogSpace.get(4)) + Parameters.minimumPossibleDecimalValue));
            }

            if(Parameters.calculateLinearFromQuadratic) {
                pn.calculateLinear();
            }
        } else if(index == 2) {
            mutateLogSpace(pn, 5, 0);
            pn.realisticBugParameters.set(index, Math.exp(pn.tauLogSpace.get(5)));

            if(Parameters.calculateLinearFromQuadratic) {
                pn.calculateLinear();
            }
        } else if(index == 3) {
            mutateLogSpace(pn, 6, Parameters.minimumAllowedMemoryLength);
            pn.realisticBugParameters.set(index, Math.exp(pn.tauLogSpace.get(6)));

            if(Parameters.calculateBugStepSize) {
                pn.calculateBugStepSize();
            }
        }
    }

    private void mutateLogSpace(ProteinNetwork pn, int index, double minimum) {
        double random;

        do {
            random = Modularity.random.nextDouble() * 0.4 - 0.2;
        } while(Math.exp(pn.tauLogSpace.get(index) + random) < minimum);

        pn.tauLogSpace.set(index, pn.tauLogSpace.get(index) + random);
    }

    private void mutateLogSpaceAllowNegative(ProteinNetwork pn, int index) {
        double random, newParameter;
        boolean conditionsAreSatisfied;

        do {
            random = Modularity.random.nextDouble() * 0.4 - 0.2;
            newParameter = Math.signum(pn.tauLogSpace.get(index) + random) * Math.exp(Math.abs(pn.tauLogSpace.get(index) + random) + Parameters.minimumPossibleDecimalValue);
            conditionsAreSatisfied = true;

            if(Math.abs(newParameter) > Parameters.sensitivity
                    || (index == 3 && newParameter < Parameters.minimumConstantAllowed)
                    || (index == 3 && newParameter > Parameters.maximumConstantAllowed)
                    || (index == 4 && newParameter < Parameters.minimumLinearAllowed)
                    || (index == 4 && newParameter > Parameters.maximumLinearAllowed)) {
                conditionsAreSatisfied = false;
            }
        } while(!conditionsAreSatisfied);

        pn.tauLogSpace.set(index, pn.tauLogSpace.get(index) + random);
    }

    private void mutateEpsilon(ProteinNetwork pn) {
        pn.epsilon += Modularity.random.nextGaussian() / 100;
        pn.calculateABFromEpsilonSensitivity();
    }

    public void fillNetwork(ProteinNetwork proteinNetwork) {
        int edgeCount = Parameters.initialProteinCount * Parameters.initialProteinCount + Parameters.initialProteinCount;
        double interactionMutationRate = Parameters.initialProteinCount * (Parameters.initialProteinCount - 1) * Parameters.mutationRateOriginal[0];
        double selfMutationRate = 2 * Parameters.initialProteinCount * Parameters.mutationRateOriginal[1];
        double sum = interactionMutationRate + selfMutationRate;
        interactionMutationRate /= sum;
        selfMutationRate /= sum;
        double mutation;
        boolean mutagenesisSuccessful = false;

        for(int i = 0; i < edgeCount; i++) {
            while(!mutagenesisSuccessful) {
                mutation = Modularity.random.nextDouble();

                if(mutation >= 0 && mutation < interactionMutationRate) {
                    if(mutateRandomInteraction(proteinNetwork)) mutagenesisSuccessful = true;
                }
                else if(mutation >= interactionMutationRate && mutation <= interactionMutationRate + selfMutationRate) {
                    if(mutateRandomSelfInteraction(proteinNetwork)) mutagenesisSuccessful = true;
                }
            }

            mutagenesisSuccessful = false;
        }
    }
}
