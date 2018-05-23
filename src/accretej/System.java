package accretej;

import javax.rmi.CORBA.Util;

public class System {
    Star centralBody;
    Planet
        planetHead = null,
        failedPlanets = null,
        escapedMoons = null;

    static final double
        B = 1.2E-5, // used in critical mass calculations
        PROTOPLANET_MASS	= 1.0E-15, // TODO: This is in stellar masses, check the validity of this. Vesta is 1.3028E-10 by comparison.
        DUST_DENSITY_COEFF = 2.0E-3,	// TODO: Read Dole's paper and figure out what "A" is.
        ALPHA = 5.0,	// TODO: Used in density calcs (how and why?)
        N = 3.0,	// TODO: Used in density calcs (how and why?)
        K = 50.0; // K = gas/dust ratio TODO: is this accurate for a protoplanetary disk?

    // planetary formation variables, these are mostly temporary
    private boolean
        dustLeft = false,
        doMoons,
        verbose,
        extraVerbose;
    public boolean habitable = false;
    private double
        radiusInner = 0.0,
        radiusOuter = 0.0,
        reducedMass = 0.0,
        dustDensity = 0.0,
        cloudEccentricity = 0.0,
        outerPlanetLimit = 0.0;
    private DustRecord dustHead = null;

    System(boolean doMoons, boolean verbose, boolean extraVerbose) {
        this.doMoons = doMoons;
        this.verbose = verbose;
        this.extraVerbose = extraVerbose;
        // centralBody = Utils.instance().randomStar();
        centralBody = Utils.instance().randomGStar();
        centralBody.setAge();
        distributePlanetaryMasses();
        migratePlanets();
        setEnvironments();
    }

    private void distributePlanetaryMasses() {
        double sma, ecc, mass, criticalMass, innerBound, outerBound;
        WDouble dustMass = new WDouble(0.0), gasMass = new WDouble(0.0);

        setInitialConditions();
        innerBound = centralBody.innermostPlanet();

        if(outerPlanetLimit == 0.0) {
            outerBound = centralBody.outermostPlanet();
        } else {
            outerBound = outerPlanetLimit;
        }
        
        while(dustLeft) {
            sma = Utils.instance().randomNumber(innerBound, outerBound);
            ecc = Utils.instance().randomEccentricity();

            mass = PROTOPLANET_MASS;

            if(verbose && extraVerbose) {
                // this really isn't significant
                java.lang.System.out.println("Checking at " + String.format("%1$,.2f", sma) + " AU...");
            }

            if(dustAvailable(innerEffectLimit(sma, ecc, mass), outerEffectLimit(sma, ecc, mass))) {
                if(verbose) {
                    java.lang.System.out.println("Injecting protoplanet at " + String.format("%1$,.2f", sma) + " AU...");
                }
                dustDensity = DUST_DENSITY_COEFF * Math.sqrt(centralBody.mass) * Math.exp(-ALPHA * Math.pow(sma, (1.0 / N)));
                criticalMass = criticalMass(sma, ecc);
                mass = accreteDust(mass, dustMass, gasMass, sma, ecc, criticalMass, innerBound, outerBound);

                dustMass.value += PROTOPLANET_MASS;

                if(mass > PROTOPLANET_MASS) {
                    coalescePlanetesimals(sma, ecc, mass, criticalMass, dustMass.value, gasMass.value, innerBound, outerBound);
                } else {
                    // save the failed planetesimal
                    Planet p = new Planet();
                    p.sma = sma;
                    p.eccentricity = ecc;
                    p.mass = dustMass.value + gasMass.value;
                    p.dustMass = dustMass.value;
                    p.gasMass = gasMass.value;
                    if(failedPlanets == null) {
                        failedPlanets = p;
                    } else {
                        failedPlanets.append(p);
                    }

                    if(verbose) {
                        java.lang.System.out.println("Failed due to larger neighbor.");
                    }
                }
            } else if(verbose && extraVerbose) {
                // this really isn't significant.
                java.lang.System.out.println("Failed due to lack of dust.");
            }
        }
    }

    private double innerEffectLimit(double sma, double ecc, double mass) {
        return (sma * (1.0 - ecc) * (1.0 - mass) / (1.0 + cloudEccentricity)); // TODO: Figure out what this means to the process and adjust as needed.
    }

    private double outerEffectLimit(double sma, double ecc, double mass) {
        return (sma * (1.0 + ecc) * (1.0 + mass) / (1.0 - cloudEccentricity)); // TODO: Figure out what this means to the process and adjust as needed.
    }

    private double criticalMass(double sma, double ecc) {
        double periapsis = SystemObject.periapsis(sma, ecc);
        double temp = periapsis * Math.sqrt(centralBody.luminosity);
        return B * Math.pow(temp, -0.75); // TODO: why is B the value it is?
    }

    private boolean dustAvailable(double innerRange, double outerRange) {
        DustRecord currentBand;
        boolean dustHere = false;
        currentBand = this.dustHead;

        while(currentBand != null && currentBand.outerEdge < innerRange) {
            currentBand = currentBand.next;
        }
        if(currentBand != null) {
            dustHere = currentBand.dustPresent;
        }
        while(currentBand != null && currentBand.innerEdge < outerRange) {
            dustHere = dustHere || currentBand.dustPresent;
            currentBand = currentBand.next;
        }

        return dustHere;
    }

    private double accreteDust(double seedMass, WDouble newDust, WDouble newGas, double sma, double ecc, double criticalMass, double innerBound, double outerBound) {
        double newMass = seedMass, tempMass;

        do {
            tempMass = newMass;
            newMass = collectDust(newMass, newDust, newGas, sma, ecc, criticalMass, dustHead);
        } while(!(newMass - tempMass < 0.0001 * tempMass));

        seedMass += newMass;
        updateDustLanes(radiusInner, radiusOuter, seedMass, criticalMass, innerBound, outerBound);

        return seedMass;
    }

    private double collectDust(double lastMass, WDouble newDust, WDouble newGas, double sma, double ecc, double criticalMass, DustRecord band) {
        double massDensity, temp, temp1, temp2, tempDensity, bandWidth, width, volume, gasDensity = 0.0, newMass, nextMass;
        WDouble nextDust = new WDouble(0.0), nextGas = new WDouble(0.0);

        temp = lastMass / (1.0 + lastMass);
        reducedMass = Math.pow(temp, (1.0 / 4.0));
        radiusInner = innerEffectLimit(sma, ecc, reducedMass);
        radiusOuter = outerEffectLimit(sma, ecc, reducedMass);

        if(radiusInner < 0.0) {
            radiusInner = 0.0;
        }

        if(band == null) {
            return 0.0;
        } else {
            if(!band.dustPresent) {
                tempDensity = 0.0;
            } else {
                tempDensity = dustDensity;
            }

            if(lastMass < criticalMass || !band.gasPresent) {
                massDensity = tempDensity;
            } else {
                massDensity = K * tempDensity / (1.0 + Math.sqrt(criticalMass / lastMass) * (K - 1.0));
                gasDensity = massDensity - tempDensity;
            }

            if(band.outerEdge <= radiusInner || band.innerEdge >= radiusOuter) {
                return collectDust(lastMass, new WDouble(newDust.value), new WDouble(newGas.value), sma, ecc, criticalMass, band.next);
            } else {
                bandWidth = (radiusOuter - radiusInner);

                temp1 = radiusOuter - band.outerEdge;
                if (temp1 < 0.0)
                    temp1 = 0.0;
                width = bandWidth - temp1;

                temp2 = band.innerEdge - radiusInner;
                if (temp2 < 0.0)
                    temp2 = 0.0;
                width = width - temp2;

                temp = 4.0 * Math.PI * Math.pow(sma,2.0) * reducedMass
                    * (1.0 - ecc * (temp1 - temp2) / bandWidth);
                volume = temp * width;

                newMass  = volume * massDensity;
                newGas.value  = volume * gasDensity;
                newDust.value = newMass - newGas.value;

                nextMass = collectDust(lastMass, nextDust, nextGas, sma, ecc, criticalMass, band.next);

                newGas.value += nextGas.value;
                newDust.value += nextDust.value;

                return newMass + nextMass;
            }
        }
    }

    private void updateDustLanes(double min, double max, double mass, double criticalMass, double innerBound, double outerBound) {
        boolean gas;
        DustRecord node1, node2, node3;

        dustLeft = false;
        if(mass > criticalMass) {
            gas = false;
        } else {
            gas = true;
        }
        node1 = dustHead;
        while(node1 != null) {
            if(node1.innerEdge < min && node1.outerEdge > max) {
                node2 = new DustRecord();
                node2.innerEdge = min;
                node2.outerEdge = max;
                if(node1.gasPresent) {
                    node2.gasPresent = gas;
                } else {
                    node2.gasPresent = false;
                }
                node2.dustPresent = false;
                node3 = new DustRecord();
                node3.innerEdge = max;
                node3.outerEdge = node1.outerEdge;
                node3.gasPresent = node1.gasPresent;
                node3.dustPresent = node1.dustPresent;
                node3.next = node1.next;
                node1.next = node2;
                node2.next = node3;
                node1.outerEdge = min;
                node1 = node3.next;
            } else {
                if(node1.innerEdge < max && node1.outerEdge > max) {
                    node2 = new DustRecord();
                    node2.next = node1.next;
                    node2.dustPresent = node1.dustPresent;
                    node2.gasPresent = node1.gasPresent;
                    node2.outerEdge = node1.outerEdge;
                    node2.innerEdge = max;
                    node1.next = node2;
                    node1.outerEdge = max;
                    if(node1.gasPresent) {
                        node1.gasPresent = gas;
                    } else {
                        node1.gasPresent = false;
                    }
                    node1.dustPresent = false;
                    node1 = node2.next;
                } else {
                    if(node1.innerEdge < min && node1.outerEdge > min) {
                        node2 = new DustRecord();
                        node2.next = node1.next;
                        node2.dustPresent = false;
                        if(node1.gasPresent) {
                            node2.gasPresent = gas;
                        } else {
                            node2.gasPresent = false;
                        }
                        node2.outerEdge = node1.outerEdge;
                        node2.innerEdge = min;
                        node1.next = node2;
                        node1.outerEdge = min;
                        node1 = node2.next;
                    } else {
                        if(node1.innerEdge >= min && node1.outerEdge <= max) {
                            if(node1.gasPresent) {
                                node1.gasPresent = gas;
                            }
                            node1.dustPresent = false;
                            node1 = node1.next;
                        } else {
                            if(node1.outerEdge < min || node1.innerEdge > max) {
                                node1 = node1.next;
                            }
                        }
                    }
                }
            }
        }
        node1 = dustHead;
        while(node1 != null) {
            if(node1.dustPresent && node1.outerEdge >= innerBound && node1.innerEdge <= outerBound) {
                dustLeft = true;
            }
            node2 = node1.next;
            if(node2 != null) {
                if(node1.dustPresent == node2.dustPresent && node1.gasPresent == node2.gasPresent) {
                    node1.outerEdge = node2.outerEdge;
                    node1.next = node2.next;
                }
            }
            node1 = node1.next;
        }
    }

    public void coalescePlanetesimals(double sma, double ecc, double mass, double criticalMass, double dustMass, double gasMass, double innerBound, double outerBound) {
        Planet thePlanet = null, nextPlanet = null, prevPlanet = null;
        boolean finished = false;
        double temp, diff, dist1, dist2;

// First we try to find an existing planet with an over-lapping orbit.
        for(thePlanet = planetHead; thePlanet != null; thePlanet = thePlanet.next) {
            diff = thePlanet.sma - sma;

            if(diff > 0.0) {
                dist1 = (sma * (1.0 + ecc) * (1.0 + reducedMass)) - sma;
                // apoapsis
                reducedMass = Math.pow(thePlanet.mass / (1.0 + thePlanet.mass), 1.0 / 4.0);
                dist2 = thePlanet.sma - (thePlanet.sma * (1.0 - thePlanet.eccentricity) * (1.0 - reducedMass));
            } else {
                dist1 = sma - (sma * (1.0 - ecc) * (1.0 - reducedMass));
                // periapsis
                reducedMass = Math.pow(thePlanet.mass / (1.0 + thePlanet.mass), 1.0 / 4.0);
                dist2 = (thePlanet.sma * (1.0 + thePlanet.eccentricity) * (1.0 + reducedMass)) - thePlanet.sma;
            }

            if(Math.abs(diff) <= Math.abs(dist1) || Math.abs(diff) <= Math.abs(dist2)) {
                WDouble newDust = new WDouble(0.0), newGas = new WDouble(0.0);
                double newSMA = (thePlanet.mass + mass) / ((thePlanet.mass / thePlanet.sma) + (mass / sma));

                temp = thePlanet.mass * Math.sqrt(thePlanet.sma) * Math.sqrt(1.0 - Math.pow(thePlanet.eccentricity, 2.0));
                temp = temp + (mass * Math.sqrt(sma) * Math.sqrt(Math.sqrt(1.0 - Math.pow(ecc, 2.0))));
                temp = temp / ((thePlanet.mass + mass) * Math.sqrt(newSMA));
                temp = 1.0 - Math.pow(temp, 2.0);
                if(temp < 0.0 || temp >= 1.0) {
                    temp = 0.0;
                }
                ecc = Math.sqrt(temp);

                if(doMoons) {
                    double existingMass = 0.0;

                    if(thePlanet.moonHead != null) {
                        Planet m;

                        for(m = thePlanet.moonHead; m != null; m = m.next) {
                            existingMass += m.mass;
                        }
                    }

                    if(mass < criticalMass) {
                        if(SystemObject.massInEarthMasses(mass) < 2.5 && SystemObject.massInEarthMasses(mass) > .0001 && existingMass < thePlanet.mass * .05) {
                            Planet theMoon = new Planet();

                            theMoon.sma = sma;
                            theMoon.eccentricity = ecc;
                            theMoon.mass = mass;
                            theMoon.dustMass = dustMass;
                            theMoon.gasMass = gasMass;
                            theMoon.isMoon = true;

                            if(theMoon.dustMass + theMoon.gasMass > thePlanet.dustMass + thePlanet.gasMass) {
                                double tempDust = thePlanet.dustMass;
                                double tempGas  = thePlanet.gasMass;
                                double tempMass = thePlanet.mass;

                                thePlanet.dustMass = theMoon.dustMass;
                                thePlanet.gasMass = theMoon.gasMass;
                                thePlanet.mass = theMoon.mass;

                                theMoon.dustMass = tempDust;
                                theMoon.gasMass = tempGas;
                                theMoon.mass = tempMass;
                            }

                            if(thePlanet.moonHead == null) {
                                thePlanet.moonHead = theMoon;
                            } else {
                                theMoon.next = thePlanet.moonHead;
                                thePlanet.moonHead = theMoon;
                            }

                            finished = true;

                            if(verbose) {
                                java.lang.System.out.println("Moon Captured... " + String.format("%1$,.2f", thePlanet.sma) + "AU (" +
                                    String.format("%1$,.2f", thePlanet.massInEarthMasses()) + "EM) <- " + String.format("%1$,.2f", Planet.massInEarthMasses(mass)) + "EM");
                            }
                        } else {
                            // save the moon
                            Planet p = new Planet();
                            p.sma = sma;
                            p.eccentricity = ecc;
                            p.mass = dustMass + gasMass;
                            p.dustMass = dustMass;
                            p.gasMass = gasMass;
                            if(escapedMoons == null) {
                                escapedMoons = p;
                            } else {
                                escapedMoons.append(p);
                            }

                            if(verbose) {
                                String msg = "Moon Escapes... " + String.format("%1$,.2f", thePlanet.sma) + " AU (" + String.format("%1$,.2f", thePlanet.massInEarthMasses()) + "EM)";
                                if(existingMass >= thePlanet.mass * 0.05) {
                                    msg += " (big moons)";
                                }
                                if(Planet.massInEarthMasses(mass) >= 2.5) {
                                    msg += ", too big";
                                } else if(Planet.massInEarthMasses(mass) <= 0.0001) {
                                    msg += ", too small";
                                }
                                msg += " moon at " + String.format("%1$,.2f", Planet.massInEarthMasses(mass)) + "EM";
                                java.lang.System.out.println(msg);
                            }
                        }
                    }
                }

                if(!finished) {
                    if(verbose) {
                        java.lang.System.out.println("Collision between two planetesimals: " + String.format("%1$,.2f", thePlanet.sma) + " AU (" +
                            String.format("%1$,.2f", thePlanet.massInEarthMasses()) + "EM), " + String.format("%1$,.2f", sma) + " AU (" +
                            String.format("%1$,.2f", SystemObject.massInEarthMasses(mass)) + "EM = " + String.format("%1$,.2f", SystemObject.massInEarthMasses(dustMass)) +
                            "EM dust + " + String.format("%1$,.2f", SystemObject.massInEarthMasses(gasMass)) + "EM gas [" +
                            String.format("%1$,.2f", SystemObject.massInEarthMasses(criticalMass)) + "EM])-> " + String.format("%1$,.2f", newSMA) + " AU (" + String.format("%1$,.2f", ecc) + ")");
                    }
                    temp = thePlanet.mass + mass;
                    temp = accreteDust(temp, newDust, newGas, newSMA, ecc, centralBody.luminosity, innerBound, outerBound);

                    thePlanet.sma = newSMA;
                    thePlanet.eccentricity = ecc;
                    thePlanet.mass = temp;
                    thePlanet.dustMass += dustMass + newDust.value;
                    thePlanet.gasMass += gasMass + newGas.value;
                    if(thePlanet.mass >= criticalMass) {
                        thePlanet.gasGiant = true;
                    }
                    while(thePlanet.next != null && thePlanet.next.sma < newSMA) {
                        nextPlanet = thePlanet.next;

                        if(thePlanet == planetHead) {
                            planetHead = nextPlanet;
                        } else {
                            prevPlanet.next = nextPlanet;
                        }

                        thePlanet.next = nextPlanet.next;
                        nextPlanet.next = thePlanet;
                        prevPlanet = nextPlanet;
                    }
                }

                finished = true;
                break;
            } else {
                prevPlanet = thePlanet;
            }
        }

        if(!finished) { // Planetesimals didn't collide. Make it a planet.
            thePlanet = new Planet();

            thePlanet.sma = sma;
            thePlanet.eccentricity = ecc;
            thePlanet.mass = mass;
            thePlanet.dustMass = dustMass;
            thePlanet.gasMass = gasMass;
            if(mass >= criticalMass) {
                thePlanet.gasGiant = true;
            }

            if(planetHead == null) {
                planetHead = thePlanet;
                thePlanet.next = null;
            } else if(sma < planetHead.sma) {
                thePlanet.next = planetHead;
                planetHead = thePlanet;
            } else if(planetHead.next == null) {
                planetHead.next = thePlanet;
                thePlanet.next = null;
            } else {
                nextPlanet = planetHead;
                while(nextPlanet != null && nextPlanet.sma < sma) {
                    prevPlanet = nextPlanet;
                    nextPlanet = nextPlanet.next;
                }
                thePlanet.next = nextPlanet;
                prevPlanet.next = thePlanet;
            }
        }
    }

    private void migratePlanets() {

    }

    private void setEnvironments() {
        Planet p = planetHead;
        while(p != null) {
            p.finalize(this.centralBody, doMoons);
            p = p.next;
        }

        p = planetHead;
        while(p != null) {
            if(p.habitable || p.habitableMoon) {
                this.habitable = true;
                return;
            }
            p = p.next;
        }
    }

    private void setInitialConditions() {
        this.dustHead = new DustRecord();
        this.dustHead.dustPresent = true;
        this.dustHead.gasPresent = true;
        this.dustHead.outerEdge = centralBody.stellarDustLimit();
        this.dustLeft = true;
        this.cloudEccentricity = 0.2; // TODO: not entirely sure what this represents in the process, come back to it.
    }

    public int numberOfFailedPlanets() {
        return Utils.instance().countPlanets(failedPlanets);
    }

    public int numberOfEscapedMoons() {
        return Utils.instance().countPlanets(escapedMoons);
    }

    public double massOfFailedPlanets() {
        return Utils.instance().sumMassOfPlanets(failedPlanets);
    }

    public double massOfEscapedMoons() {
        return Utils.instance().sumMassOfPlanets(escapedMoons);
    }

    public String toString() {
        String str = "New System\n";
        str = str.concat("Primary: " + centralBody.toString() + "\n");
        str = str.concat("Failed planets: " + numberOfFailedPlanets() + " (" + String.format("%1$,.4f", SystemObject.massInEarthMasses(massOfFailedPlanets())) + "EM)\n");
        str = str.concat("Escaped Moons: " + numberOfEscapedMoons() + " (" + String.format("%1$,.4f", SystemObject.massInEarthMasses(massOfEscapedMoons())) + "EM)\n");
        Planet temp = planetHead;
        int i = 1;
        while(temp != null) {
            str = str.concat("  Planet " + i++ + ": " + temp + "\n");
            temp = temp.next;
        }
        return str;
    }

    private class DustRecord {
        DustRecord next = null;
        double innerEdge = 0.0, outerEdge = 0.0;
        boolean dustPresent = false, gasPresent = false;
    }

    private class WDouble {
        double value = 0.0;
        WDouble(double d) {
            value = d;
        }
    }
}
