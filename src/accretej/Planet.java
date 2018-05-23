package accretej;

import java.util.Iterator;
import java.util.Vector;

public class Planet extends SystemObject {
    protected double
        axialTilt = 0.0,
        radius = 0.0, // in km
        dustMass = 0.0,
        gasMass = 0.0,
        density = 0.0, // in g/cc
        orbitalPeriod = 0.0,
        coreRadius = 0.0,
        dayLength = 0.0,
        exosphericTemperature = 0.0,
        estimatedTemperature = 0.0,
        estimatedTerrestrialTemperature = 0.0,
        rmsVelocity = 0.0,
        escapeVelocity = 0.0,
        volatileGasInventory = 0.0,
        surfaceAcceleration = 0.0,
        surfaceGravity = 0.0,
        surfacePressure = 0.0,
        surfaceTemperature = 0.0,
        highTemperature = 0.0,
        lowTemperature = 0.0,
        maxTemperature = 0.0,
        minTemperature = 0.0,
        boilingPoint = 0.0,
        greenhouseRise = 0.0,
        minimumMolecularWeight = 0.0,
        hydrosphere = 0.0,
        cloudCover = 0.0,
        albedo = 0.0,
        iceCover = 0.0;
    protected int orbitalZone = -1;
    protected boolean
        gasGiant = false,
        habitableJovian = false,
        isMoon = false,
        resonantPeriod = false,
        greenhouseEffect = false,
        habitable = false,
        earthlike = false,
        habitableMoon = false;
    protected planetType type;
    protected Planet
        next = null,
        moonHead = null;
    protected Vector<AtmosphericChemical> atmosphere = new Vector<AtmosphericChemical>();

    protected enum atmosphereType {
        NONE,
        BREATHABLE,
        UNBREATHABLE,
        POISONOUS
    }

    protected enum planetType {
        tUnknown,
        tSubSubGasGiant,
        tSubGasGiant,
        tGasGiant,
        tRock,
        tVenusian,
        tTerrestrial,
        tMartian,
        tWater,
        tIce,
        tAsteroids,
        tTidallyLocked
    }

    /**
     * Performs finalizing calculation on the planet based on its given parameters.
     * Make all necessary changes to sma, eccentricity, various masses, etc... before running this.
     */
    public void finalize(Star primary, boolean doMoons) {
        this.surfaceTemperature = 0.0;
        this.highTemperature = 0.0;
        this.lowTemperature = 0.0;
        this.maxTemperature = 0.0;
        this.minTemperature = 0.0;
        this.greenhouseRise = 0.0;
        this.resonantPeriod = false;
        this.orbitalZone = orbitalZone(primary.luminosity);
        this.orbitalPeriod = orbitalPeriod(primary);
        this.axialTilt = axialTilt();
        this.exosphericTemperature = EARTH_EXOSPHERIC_TEMPERATURE / Math.pow(this.sma / primary.radiusEcosphere, 2.0);
        this.rmsVelocity = rmsVelocity(MOLECULAR_NITROGEN);
        this.coreRadius = kothariRadius();

        // Calculate the radius as a gas giant, to verify it will retain gas. Then if mass > Earth,
        // it's at least 5% gas, and retains He, it's some flavor of gas giant.

        this.density = empiricalDensity(primary.radiusEcosphere);
        this.radius = volumeRadius();
        this.surfaceAcceleration = gravitationalAcceleration();
        this.surfaceGravity = gravity();
        this.minimumMolecularWeight = minimumMolecularWeight(primary);

        if(massInEarthMasses() > 1.0 && this.gasMass / this.mass > 0.05 && this.minimumMolecularWeight <= 4.0) {
            // Gas Giant Planet
            if(this.gasMass / this.mass < 0.20) {
                this.type = planetType.tSubSubGasGiant;
            } else if(massInEarthMasses() < 20.0) {
                this.type = planetType.tSubGasGiant;
            } else {
                this.type = planetType.tGasGiant;
            }
        } else {
            // Rocky Planet
            this.radius = kothariRadius();
            this.density = volumeDensity();
            this.surfaceAcceleration = gravitationalAcceleration();
            this.surfaceGravity = gravity();

            if(this.gasMass / this.mass > 0.000001) {
                double h2Mass = this.gasMass * 0.85;
                double heMass = (this.gasMass - h2Mass) * 0.999;
                double h2Life = gasLife(MOLECULAR_HYDROGEN);
                double heLife = gasLife(HELIUM);
                double h2Loss, heLoss;

                if(h2Life < primary.age) {
                    h2Loss = (1.0 - (1.0 / Math.exp(primary.age / h2Life))) * h2Mass;
                    this.gasMass -= h2Loss;
                    this.mass -= h2Loss;
                    this.surfaceAcceleration = gravitationalAcceleration();
                    this.surfaceGravity = gravity();
                }

                if(heLife < primary.age) {
                    heLoss = (1.0 - (1.0 / Math.exp(primary.age / heLife))) * heMass;
                    this.gasMass -= heLoss;
                    this.mass -= heLoss;
                    this.surfaceAcceleration = gravitationalAcceleration();
                    this.surfaceGravity = gravity();
                }

                /* TODO: Logging singleton
                if (((h2Loss + heLoss) > .000001) && (flag_verbose & 0x0080))
                    fprintf (stderr, "%s\tLosing gas: H2: %5.3Lf EM, He: %5.3Lf EM\n", planet_id, h2Loss * SUN_MASS_IN_EARTH_MASSES, heLoss * SUN_MASS_IN_EARTH_MASSES);
                */
            }
        }

        this.dayLength = dayLength(primary);
        this.escapeVelocity = escapeVelocity();

        if(this.type == planetType.tGasGiant || this.type == planetType.tSubGasGiant || this.type == planetType.tSubSubGasGiant) {
            this.greenhouseEffect = false;
            this.volatileGasInventory = Double.MAX_VALUE;
            this.surfacePressure = Double.MAX_VALUE;
            this.boilingPoint = Double.MAX_VALUE;
            this.surfaceTemperature = Double.MAX_VALUE;
            this.greenhouseRise = 0.0;
            this.albedo = Utils.instance().about(GAS_GIANT_ALBEDO,0.1);
            this.hydrosphere = 1.0;
            this.cloudCover = 1.0;
            this.iceCover = 0.0;
            this.surfaceGravity = gravity();
            this.minimumMolecularWeight = minimumMolecularWeight(primary);
            this.surfaceGravity = Double.MAX_VALUE;
            this.estimatedTemperature = estimatedTemperature(primary);
            this.estimatedTerrestrialTemperature = estimatedTerrestrialTemperature(primary);

            if(this.estimatedTerrestrialTemperature >= FREEZING_POINT_OF_WATER && this.estimatedTerrestrialTemperature <= EARTH_AVERAGE_KELVIN + 10.0 && primary.age > 2.0E9) {
                this.habitableJovian = true;

                /* TODO: Logging singleton
                if(flag_verbose & 0x8000) {
                    fprintf (stderr, "%s\t%s (%4.2LfEM %5.3Lf By)%s with earth-like temperature (%.1Lf C, %.1Lf F, %+.1Lf C Earth).\n",
                        planet_id,
                        planet->type == tGasGiant ? "Jovian" :
                            planet->type == tSubGasGiant ? "Sub-Jovian" :
                                planet->type == tSubSubGasGiant ? "Gas Dwarf" :
                                    "Big",
                        planet->mass * SUN_MASS_IN_EARTH_MASSES,
                        primary.age /1.0E9,
                        planet->first_moon == NULL ? "" : " WITH MOON",
                        temp - FREEZING_POINT_OF_WATER,
                        32 + ((temp - FREEZING_POINT_OF_WATER) * 1.8),
                        temp - EARTH_AVERAGE_KELVIN);
                }
                */
            }
        } else {
            this.estimatedTemperature = estimatedTemperature(primary);
            this.estimatedTerrestrialTemperature = estimatedTerrestrialTemperature(primary);
            this.surfaceGravity = gravity();
            this.minimumMolecularWeight = minimumMolecularWeight(primary);
            this.greenhouseEffect = greenhouse(primary);
            this.volatileGasInventory= volatileGasInventory(primary);
            this.surfacePressure = pressure();

            if(this.surfacePressure <= 0.00001) { // TODO: May have to come back to this, was originally 0.0
                this.boilingPoint = 0.0;
            } else {
                this.boilingPoint = boilingPoint();
            }

            iterateSurfaceTemperature(primary);

            if(this.maxTemperature >= FREEZING_POINT_OF_WATER && this.minTemperature <= this.boilingPoint) {
                calculateGases(primary);
            }

            if(this.surfacePressure < 1.0) {
                if(!this.isMoon && massInEarthMasses() < ASTEROID_MASS_LIMIT) {
                    this.type = planetType.tAsteroids;
                } else {
                    this.type = planetType.tRock;
                }
            } else if(this.surfacePressure > 6000.0 && this.minimumMolecularWeight <= 2.0) { // Retains Hydrogen
                this.type = planetType.tSubSubGasGiant;
                this.atmosphere = null;
            } else { // Atmospheres:
                if(secondsToHoursRounded(this.dayLength) == secondsToHoursRounded(this.orbitalPeriod) || this.resonantPeriod) {
                    this.type = planetType.tTidallyLocked;
                } else if(this.hydrosphere >= 0.95) {
                    this.type = planetType.tWater; // >95% water
                } else if(this.iceCover >= 0.95) {
                    this.type = planetType.tIce; // >95% ice
                } else if(this.hydrosphere > 0.05) {
                    this.type = planetType.tTerrestrial; // Terrestrial
                    // else <5% water
                } else if(this.maxTemperature > this.boilingPoint) {
                    this.type = planetType.tVenusian; // Hot = Venusian
                } else if(this.gasMass / this.mass > 0.0001) { // Accreted gas
                    this.type = planetType.tIce; // But no Greenhouse
                    this.iceCover = 1.0; // or liquid water
                } else if(surfacePressure <= 250.0) { // Thin air = Martian
                    this.type = planetType.tMartian;
                } else if(this.surfaceTemperature < FREEZING_POINT_OF_WATER) {
                    this.type = planetType.tIce;
                } else {
                    this.type = planetType.tUnknown;
                    /* TODO: Logging singleton
                    if (flag_verbose & 0x0001)
                        fprintf (stderr, "%12s\tp=%4.2Lf\tm=%4.2Lf\tg=%4.2Lf\tt=%+.1Lf\t%s\t Unknown %s\n",
                            type_string (this.type),
                            this.surf_pressure,
                            this.mass * SUN_MASS_IN_EARTH_MASSES,
                            this.surf_grav,
                            this.surf_temp  - EARTH_AVERAGE_KELVIN,
                            planet_id,
                            ((int)this.day == (int)(this.orb_period * 24.0) ||
                                (this.resonant_period)) ? "(1-Face)" : ""
                        );
                    */
                }
            }
        }

        if(doMoons && !this.isMoon) {
            if(this.moonHead != null) {
                Planet p = this.moonHead;

                while(p != null) {
                    if(p.massInEarthMasses() > .000001) {
                        p.sma = this.sma;
                        p.eccentricity = this.eccentricity;
                        p.finalize(primary, false);

                        double rocheLimit = 2.44 * this.radius * Math.pow(this.density / p.density, 1.0 / 3.0);
                        double hillSphere = this.sma * KM_PER_AU * Math.pow(this.mass / (3.0 * primary.mass), 1.0 / 3.0);

                        if(rocheLimit * 3.0 < hillSphere) {
                            p.asMoonSMA = Utils.instance().randomNumber(rocheLimit * 1.5, hillSphere / 2.0) / KM_PER_AU;
                            p.asMoonEccentricty = Utils.instance().randomEccentricity();
                        }

                        /* TODO: logging singleton
                        if(flag_verbose & 0x40000) {
                            fprintf (stderr,
                                "   Roche limit: R = %4.2Lg, rM = %4.2Lg, rm = %4.2Lg -> %.0Lf km\n"
                                "   Hill Sphere: a = %4.2Lg, m = %4.2Lg, M = %4.2Lg -> %.0Lf km\n"
                                "%s Moon orbit: a = %.0Lf km, e = %.0Lg\n",
                                this.radius, this.density, ptr->density,
                                roche_limit,
                                this.a * KM_PER_AU, this.mass * SOLAR_MASS_IN_KILOGRAMS, sun->mass * SOLAR_MASS_IN_KILOGRAMS,
                                hill_sphere,
                                moon_id,
                                ptr->moon_a * KM_PER_AU, ptr->moon_e
                            );
                        }

                        if(flag_verbose & 0x1000) {
                            fprintf (stderr, "  %s: (%7.2LfEM) %d %4.2LgEM\n",
                                planet_id,
                                this.mass * SUN_MASS_IN_EARTH_MASSES,
                                n,
                                ptr->mass * SUN_MASS_IN_EARTH_MASSES);
                        }
                        */
                        if(p.habitable) {
                            this.habitableMoon = true;
                        }
                    }
                    p = p.next;
                }
            }
        }

        // ensure we don't have a mislabeled planet.
        if(this.type != planetType.tSubSubGasGiant && this.type != planetType.tSubGasGiant && this.type != planetType.tGasGiant) {
            this.gasGiant = false;
        } else {
            this.gasGiant = true;
            return; // no need to continue
        }

        // check for those golden planets
        if(breathableAtmosphere() == atmosphereType.BREATHABLE) {
            if(!(this.resonantPeriod || secondsToHoursRounded(this.dayLength) == secondsToHoursRounded(this.orbitalPeriod))) {
                // TODO: Investigate habitable tidally locked planets.
                this.habitable = true;

                double relTemp = this.surfaceTemperature - FREEZING_POINT_OF_WATER - EARTH_AVERAGE_CELSIUS;
                double pressure = this.surfacePressure / EARTH_SURF_PRES_IN_MILLIBARS;

                if(this.surfaceGravity >= 0.8 && this.surfaceGravity <= 1.2 && relTemp >= -2.0 && relTemp <= 3.0 && this.iceCover <= 0.1 && pressure >= 0.5 && pressure <= 2.0 && this.cloudCover >= 0.4 && this.cloudCover <= 0.8 && this.hydrosphere >= 0.5 && this.hydrosphere <= 0.8 && this.type != planetType.tWater) {
                    this.earthlike = true;
                }
            }
        }
    }

    public void calculateGases(Star primary) {
        if(this.surfacePressure > 0.0) {
            Chemical[] chems = Utils.instance().getChemicals();
            double amount[] = new double[chems.length];
            double totalAmount = 0.0;
            double pressure = this.surfacePressure / 1000.0; // convert from millibars to bars

            for(int i = 0; i < chems.length; i++) {
                double yp = chems[i].boilingPoint / (373.0 * ((Math.log(pressure + 0.001) / -5050.5) + (1.0 / 373.0)));

                if(yp >= 0 && yp < this.lowTemperature && chems[i].weight >= this.minimumMolecularWeight) {
                    double vrms = rmsVelocity(chems[i].weight);
                    double pvrms = Math.pow(1.0 / (1.0 + vrms / this.escapeVelocity), primary.age / 1e9);
                    double abund = chems[i].abundanceS;
                    double react = 1.0;
                    double fract = 1.0;
                    double pres2 = 1.0;

                    if(chems[i].atomicNumber == 18) {
                        react = .15 * primary.age/4e9;
                    } else if(chems[i].atomicNumber == 2) {
                        abund = abund * (0.001 + (this.gasMass / this.mass));
                        pres2 = (0.75 + pressure);
                        react = Math.pow(1.0 / (1.0 + chems[i].reactivity), primary.age / 2e9 * pres2);
                    } else if((chems[i].atomicNumber == 8 || chems[i].atomicNumber == 912) && primary.age > 2e9 && this.surfaceTemperature > 270.0 && this.surfaceTemperature < 400.0) {
                        /*	pres2 = (0.65 + pressure/2); Breathable - M: .55-1.4 */
                        pres2 = (0.89 + pressure / 4); /*	Breathable - M: .6 -1.8 */
                        react = Math.pow(1.0 / (1.0 + chems[i].reactivity), Math.pow(primary.age / 2e9, 0.25) * pres2);
                    } else if(chems[i].atomicNumber == 902 && primary.age > 2e9 && this.surfaceTemperature > 270 && this.surfaceTemperature < 400) {
                        pres2 = (0.75 + pressure);
                        react = Math.pow(1.0 / (1.0 + chems[i].reactivity), Math.pow(primary.age / 2e9, 0.5) * pres2);
                        react *= 1.5;
                    } else {
                        pres2 = (0.75 + pressure);
                        react = Math.pow(1.0 / (1.0 + chems[i].reactivity), primary.age / 2e9 * pres2);
                    }

                    fract = 1.0 - (this.minimumMolecularWeight / chems[i].weight);
                    amount[i] = abund * pvrms * react * fract;

/* TODO: Logging
                    if ((flag_verbose & 0x4000) &&
                        (strcmp(gases[i].symbol, "O") == 0 ||
                            strcmp(gases[i].symbol, "N") == 0 ||
                            strcmp(gases[i].symbol, "Ar") == 0 ||
                            strcmp(gases[i].symbol, "He") == 0 ||
                            strcmp(gases[i].symbol, "CO2") == 0))
                    {
                        fprintf (stderr, "%-5.2Lf %-3.3s, %-5.2Lf = a %-5.2Lf * p %-5.2Lf * r %-5.2Lf * p2 %-5.2Lf * f %-5.2Lf\t(%.3Lf%%)\n",
                            this.mass * SUN_MASS_IN_EARTH_MASSES,
                            gases[i].symbol,
                            amount[i],
                            abund,
                            pvrms,
                            react,
                            pres2,
                            fract,
                            100.0 * (this.gas_mass / this.mass)
                        );
                    }
*/

                    totalAmount += amount[i];
                } else {
                    amount[i] = 0.0;
                }
            }

            if(totalAmount > 0.0) {
                for(int i = 0, n = 0; i < chems.length; i++) {
                    if(amount[i] > 0.0) {
                        this.atmosphere.add(new AtmosphericChemical(chems[i], this.surfacePressure * amount[i] / totalAmount));

/* TODO: logging
                        if (flag_verbose & 0x2000) {
                            if ((this.atmosphere[n].num == AN_O) &&
                                inspired_partial_pressure (this.surfacePressure,
                                    this.atmosphere[n].surfacePressure)
                                    > gases[i].max_ipp)
                            {
                                fprintf (stderr, "%s\t Poisoned by O2\n",
                                    planet_id);
                            }
                        }
*/

                        n++;
                    }
                }

                // TODO: Do some sorting on the resulting vector
                // qsort(this.atmosphere, this.gases, sizeof(gas), diminishing_pressure);

/* TODO: Logging
                if (flag_verbose & 0x0010)
                {
                    fprintf (stderr, "\n%s (%5.1Lf AU) gases:\n",
                        planet_id, this.a);

                    for (i = 0; i < this.gases; i++)
                    {
                        fprintf (stderr, "%3d: %6.1Lf, %11.7Lf%%\n",
                            this.atmosphere[i].num,
                            this.atmosphere[i].surfacePressure,
                            100. * (this.atmosphere[i].surfacePressure /
                                this.surfacePressure)
                        );
                    }
                }
*/

            }
        }
    }

    public double radiusInMeters() {
        return this.radius * 1000.0;
    }

    public double radiusInEarthRadii() {
        return this.radius / EARTH_RADIUS;
    }

    public double ratioRadiusToEarth() {
        return EARTH_RADIUS / this.radius;
    }

    public int numberOfMoons() {
        return Utils.instance().countPlanets(this.moonHead);
    }

    public void append(Planet p) {
        if(next == null) {
            next = p;
        } else {
            next.append(p);
        }
    }

    /**
     * @param luminosity The stellar luminosity of the primary.
     * @return The "orbital zone" of this planet.
     */
    public int orbitalZone(double luminosity) {
        // TODO: Does this have anything to do with the frost line?
        if(this.sma < 4.0 * Math.sqrt(luminosity)) {
            return 1;
        }
        if(this.sma < 15.0 * Math.sqrt(luminosity)) {
            return 2;
        }
        return 3;
    }

    /**
     * @return The radius of this planet in km
     */
    public double volumeRadius() {
        return Math.pow(3.0 * (massInGrams() / this.density) / (4.0 * Math.PI), 1.0 / 3.0) / CM_PER_KM;
    }

    /**
     * This formula is listed as eq.9 in Fogg's article, although some typos crop up in that eq.
     * See "The Internal Constitution of Planets", by Dr. D. S. Kothari, Mon. Not. of the Royal
     * Astronomical Society, vol 96 pp.833-843, 1936 for the derivation.  Specifically, this is
     * Kothari's eq.23, which appears on page 840.
     * @return The radius of the planet in kilometers.
     */
    public double kothariRadius() {
        double temp1, temp, temp2, atomicWeight, atomicNumber;
        double A1_20 = 6.485E12, A2_20 = 4.0032E-8, BETA_20 = 5.71E12;

        if(this.orbitalZone == 1) {
            if(this.gasGiant) {
                atomicWeight = 9.5;
                atomicNumber = 4.5;
            } else {
                atomicWeight = 15.0;
                atomicNumber = 8.0;
            }
        } else if(this.orbitalZone == 2) {
            if(this.gasGiant) {
                atomicWeight = 2.47;
                atomicNumber = 2.0;
            } else {
                atomicWeight = 10.0;
                atomicNumber = 5.0;
            }
        } else {
            if(this.gasGiant) {
                atomicWeight = 7.0;
                atomicNumber = 4.0;
            } else {
                atomicWeight = 10.0;
                atomicNumber = 5.0;
            }
        }

        temp1 = atomicWeight * atomicNumber;

        temp = (2.0 * BETA_20 * Math.pow(SUN_MASS_IN_GRAMS, 1.0 / 3.0)) / (A1_20 * Math.pow(temp1, 1.0 / 3.0));
        temp2 = A2_20 * Math.pow(atomicWeight, 4.0 / 3.0) * Math.pow(SUN_MASS_IN_GRAMS, 2.0 / 3.0);
        temp2 = temp2 * Math.pow(this.mass, 2.0 / 3.0);
        temp2 = temp2 / (A1_20 * Math.pow(atomicNumber, 2.0));
        temp2 = 1.0 + temp2;
        temp = temp / temp2;
        temp = (temp * Math.pow(this.mass,(1.0 / 3.0))) / CM_PER_KM;

        temp /= 1.004; // Make Earth = actual earth, called "JIMS_FUDGE" in the original code.

        return temp;
    }

    /**
     * @param radiusEcosphere SMA of the ecosphere in AU.
     * @return the empirical density of the planet in g/cc
     */
    public double empiricalDensity(double radiusEcosphere) {
        double temp;

        temp = Math.pow(massInEarthMasses(), 1.0 / 8.0);
        temp = temp * Math.sqrt(Math.sqrt(radiusEcosphere / this.sma));
        if(this.gasGiant) {
            return temp * 1.2;
        } else {
            return temp * 5.5;
        }
    }

    /**
     * @return The volume density based on assumption the planet is a perfect sphere, in g/cc
     */
    public double volumeDensity() {
        return massInGrams() / ((4.0 * Math.PI * Math.pow(this.radius * CM_PER_KM, 3.0)) / 3.0);
    }

    /**
     * @param largerBody The larger body this planet orbits around
     * @return The orbital period of this planet in seconds.
     */
    public double orbitalPeriod(SystemObject largerBody) {
        return Math.sqrt(Math.pow(smaInMeters(), 3.0) / largerBody.mu()) * Math.PI * 2.0;
    }

    /**
     * @param separation The average distance between the two bodies.
     * @param smallerBody The smaller body in this equation.
     * @param largerBody The larger body in this equation, usually the star.
     * @return The orbital period of this binary pair in seconds.
     */
    public static double binaryOrbitalPeriod(double separation, SystemObject smallerBody, SystemObject largerBody) {
        return Math.sqrt(Math.pow(smaInMeters(separation), 3.0) / (G * (smallerBody.massInKg() + largerBody.massInKg()))) * Math.PI * 2.0;
    }

    /**
     * Fogg's information for this routine came from Dole "Habitable Planets for Man", Blaisdell
     * Publishing Company, NY, 1964. From this, he came up with his eq.12, which is the equation
     * for the 'base_angular_velocity' below. He then used an equation for the change in angular
     * velocity per time (dw/dt) from P. Goldreich and S. Soter's paper "Q in the Solar System" in
     * Icarus, vol 5, pp.375-389 (1966). Using as a comparison the change in angular velocity
     * for the Earth, Fogg has come up with an approximation for our new planet (his eq.13) and
     * take that into account. This is used to find 'changeInAngularVelocity' below.
     * @return The length of this planet's rotation in seconds.
     */
    public double dayLength(Star primary) {
        final double
            J = 1.46E-19, // Used in day-length calcs (cm2/sec2 g)
            CHANGE_IN_EARTH_ANGULAR_VELOCITY = -1.3E-15; // Units of radians/sec/year

        double k2, baseAngularVelocity, changeInAngularVelocity, angularVelocity, spinResonanceFactor, dayInSeconds;
        boolean stopped = false;

        if(this.gasGiant) {
            k2 = 0.24;
        } else {
            k2 = 0.33;
        }

        baseAngularVelocity = Math.sqrt(2.0 * J * massInGrams() / (k2 * Math.pow(radius * CM_PER_KM, 2.0)));

        // This next calculation determines how much the planet's rotation is slowed by the presence of the star.
        changeInAngularVelocity = CHANGE_IN_EARTH_ANGULAR_VELOCITY * (this.density / EARTH_DENSITY) * ((this.radius * CM_PER_KM) / EARTH_RADIUS) * (1.0 / massInEarthMasses()) * Math.pow(primary.mass, 2.0) * (1.0 / Math.pow(this.sma, 6.0));
        angularVelocity = baseAngularVelocity + (changeInAngularVelocity * primary.age);
        if(angularVelocity <= 0.0) {
            stopped = true;
            dayInSeconds = Double.MAX_VALUE;
        } else {
            dayInSeconds = 2.0 * Math.PI * angularVelocity;
        }

        // tidal or resonant locking?
        if(dayInSeconds >= this.orbitalPeriod || stopped) {
            if(this.eccentricity > 0.1) {
                spinResonanceFactor = (1.0 - this.eccentricity) / (1.0 + this.eccentricity);
                this.resonantPeriod = true;
                return spinResonanceFactor * this.orbitalPeriod;
            } else {
                return this.orbitalPeriod;
            }
        }

        return dayInSeconds;
    }

    /**
     * Randomly generates an axial tilt for this planet
     * @return The axial tilt in degrees.
     */
    public double axialTilt() {
        final double EARTH_AXIAL_TILT = 23.5;
        double temp;

        temp = Math.pow(this.sma, 0.2) * Utils.instance().about(EARTH_AXIAL_TILT,0.4);
        return temp % 360.0;
    }

    /**
     * This function implements the escape velocity calculation. Note that it appears that Fogg's
     * eq.15 is incorrect.
     * @return The escape velocity in m/s.
     */
    public double escapeVelocity() {
        return Math.sqrt((2.0 * mu()) / radiusInMeters());
    }

    /**
     * This is Fogg's eq.16.  The molecular weight (usually assumed to be N2) is used as the basis
     * of the Root Mean Square (RMS) velocity of the molecule or atom.
     * @param molecularWeight
     * @return The RMS Velocity of the molecule in m/s
     */
    public double rmsVelocity(double molecularWeight) {
        return Math.sqrt((3.0 * MOLAR_GAS_CONST * this.exosphericTemperature) / molecularWeight);
    }

    /**
     * This function returns the smallest molecular weight retained by the body, which is useful
     * for determining the atmosphere composition.
     * @return
     */
    public double molecularLimit() {
        return (3.0 * MOLAR_GAS_CONST * this.exosphericTemperature) / Math.pow((escapeVelocity() / GAS_RETENTION_THRESHOLD), 2.0);
    }

    /**
     * This function calculates the surface acceleration of a planet.
     * @return The surface acceleration in m/s.
     */
    public double gravitationalAcceleration() {
        return mu() / Math.pow(radiusInMeters(), 2.0);
    }

    /**
     * @return The gravity of this planet in Earth gravities.
     */
    public double gravity() {
        return gravitationalAcceleration() / EARTH_ACCELERATION;
    }

    /**
     * This implements Fogg's eq.17.  The 'inventory' returned is unitless.
     * @param primary The primary for this planet.
     * @return The unitless "inventory".
     */
    public double volatileGasInventory(Star primary) {
        double proportionalConst, temp;

        double velocityRatio = this.escapeVelocity / this.rmsVelocity;
        if(velocityRatio >= GAS_RETENTION_THRESHOLD) {
            switch(this.orbitalZone) {
            case 1:
                proportionalConst = 140000.0; // 100 -> 140 JLB
                break;
            case 2:
                proportionalConst = 75000.0;
                break;
            case 3:
                proportionalConst = 250.0;
                break;
            default:
                proportionalConst = 0.0;
                // printf("Error: orbital zone not initialized correctly!\n");
                break;
            }
            temp = (proportionalConst * massInEarthMasses()) / primary.mass;
            /* TODO: I'm sure the original author had something in mind here...
            temp2 = Utils.instance().about(temp1, 0.2);
            temp2 = temp; // what? */
            if(this.greenhouseEffect || this.gasMass / this.mass > 0.000001) {
                return temp;
            } else {
                return temp / 140.0; // 100 -> 140 JLB
            }
        } else {
            return 0.0;
        }
    }

    /**
     * This implements Fogg's eq.18. JLB: Aparently this assumed that earth pressure = 1000mb. I've
     * added a fudge factor (EARTH_SURF_PRES_IN_MILLIBARS / 1000.) to correct for that.
     * @return Pressure in millibars.
     */
    public double pressure() {
        return this.volatileGasInventory * this.surfaceGravity * (EARTH_SURF_PRES_IN_MILLIBARS / 1000.0) / Math.pow(ratioRadiusToEarth(), 2.0);
    }

    /**
     * This function returns the boiling point of water in this planet's atmosphere. This is Fogg's
     * eq.21.
     * @return The boiling point in units of Kelvin.
     */
    public double boilingPoint() {
        return 1.0 / ((Math.log(this.surfacePressure / 1000.0) / -5050.5) + (1.0 / 373.0));
    }

    /**
     * This function is Fogg's eq.22. Given the volatile gas inventory and planetary radius of a
     * planet (in Km), this function returns the fraction of the planet covered with water. I have
     * changed the function very slightly:	 the fraction of Earth's surface covered by water is
     * 71%, not 75% as Fogg used.
     * @return The fraction of the planet covered by water.
     */
    public double waterCoverage() {
        double temp;

        temp = (0.708 * this.volatileGasInventory / 1000.0) * Math.pow(ratioRadiusToEarth(), 2.0);
        if(temp >= 1.0) {
            return 1.0;
        } else {
            return temp;
        }
    }

    /**
     * This function returns the fraction of cloud cover available. This is Fogg's eq.23. See Hart
     * in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation. This equation is Hart's eq.3. I
     * have modified it slightly using constants and relationships from Glass's book "Introduction
     * to Planetary Geology", p.46. The 'CLOUD_COVERAGE_FACTOR' is the amount of surface area on
     * Earth covered by one Kg. of cloud.
     * @return The cloud coverage fraction available.
     */
    public double cloudFraction() {
        double waterVaporKg, fraction, surfaceArea, hydroMass;
        final double
            Q2_36 = 0.0698, // 1/Kelvin
            CLOUD_COVERAGE_FACTOR = 1.839E-8; // Km2/kg

        if(this.minimumMolecularWeight > WATER_VAPOR) {
            return 0.0;
        } else {
            surfaceArea = 4.0 * Math.PI * Math.pow(this.radius, 2.0);
            hydroMass = this.hydrosphere * surfaceArea * EARTH_WATER_MASS_PER_AREA;
            waterVaporKg = (0.00000001 * hydroMass) * Math.exp(Q2_36 * (this.surfaceTemperature - EARTH_AVERAGE_KELVIN));
            fraction = CLOUD_COVERAGE_FACTOR * waterVaporKg / surfaceArea;
            if(fraction >= 1.0) {
                return 1.0;
            } else {
                return fraction;
            }
        }
    }

    /**
     * This is Fogg's eq.24. See Hart[24] in Icarus vol.33, p.28 for an explanation. I have changed
     * a constant from 70 to 90 in order to bring it more in line with the fraction of the Earth's
     * surface covered with ice, which is approximatly .016 (=1.6%).
     * @return The fraction of the planet's surface covered by ice.
     */
    public double iceFraction() {
        double temp, surfTemp;
        surfTemp = this.surfaceTemperature;

        if(surfTemp > 328.0) {
            surfTemp = 328.0;
        }
        temp = Math.pow((328.0 - surfTemp) / 90.0, 5.0);
        if(temp > 1.5 * this.hydrosphere) {
            temp = 1.5 * this.hydrosphere;
        }
        if(temp >= 1.0) {
            return 1.0;
        } else {
            return temp;
        }
    }

    /**
     * This is Fogg's eq.19.
     * @param primary The star this planet orbits.
     * @return The temperature in Kelvin.
     */
    public double effectiveTemperature(Star primary, double givenAlbedo) {
        final double EARTH_EFFECTIVE_TEMP = 250.0; // Units of degrees Kelvin (was 255)

        return Math.sqrt(primary.radiusEcosphere / this.sma) * Math.sqrt(Math.sqrt((1.0 - givenAlbedo) / (1.0 - EARTH_ALBEDO))) * EARTH_EFFECTIVE_TEMP;
    }

    /**
     * @param primary The star this planet orbits.
     * @return The temperature in Kelvin.
     */
    public double estimatedTemperature(Star primary) {
        return Math.sqrt(primary.radiusEcosphere / this.sma) * Math.sqrt(Math.sqrt((1.0 - this.albedo) / (1.0 - EARTH_ALBEDO))) * EARTH_AVERAGE_KELVIN;
    }

    /**
     * @param primary The star this planet orbits.
     * @return The temperature in Kelvin.
     */
    public double estimatedTerrestrialTemperature(Star primary) {
        return Math.sqrt(primary.radiusEcosphere / this.sma) * EARTH_AVERAGE_KELVIN;
    }

    /**
     * Old greenhouse: Note that if the orbital radius of the planet is greater than or equal to
     * R_inner, 99% of it's volatiles are assumed to have been deposited in surface reservoirs
     * (otherwise, it suffers from the greenhouse effect).
     *
     * if((sma < r_greenhouse) && (orbitalZone == 1))
     *
     * The new definition is based on the inital surface temperature and what state water is in. If
     * it's too hot, the water will never condense out of the atmosphere, rain down and form an
     * ocean. The albedo used here was chosen so that the boundary is about the same as the old
     * method Neither zone, nor r_greenhouse are used in this version - JLB
     * @param primary The star this planet orbits
     * @return Whether a greenhouse effect exists
     */
    boolean greenhouse(Star primary) {
        if(effectiveTemperature(primary, GREENHOUSE_TRIGGER_ALBEDO) > FREEZING_POINT_OF_WATER) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of Earth's Atmosphere"
     * article. The effective temperature given is in units of Kelvin, as is the rise in
     * temperature produced by the greenhouse effect, which is returned. I tuned this by changing a
     * pow(x,.25) to pow(x,.4) to match Venus - JLB
     * @return the rise in temperature produced by the greenhouse effect in Kelvin
     */
    public double greenhouseRise(double effectiveTemperature) {
        final double EARTH_CONVECTION_FACTOR = 0.43; // from Hart, eq.20

        double convection_factor = EARTH_CONVECTION_FACTOR * Math.pow(this.surfacePressure / EARTH_SURF_PRES_IN_MILLIBARS, 0.4);
        double rise = (Math.sqrt(Math.sqrt(1.0 + 0.75 * opticalDepth())) - 1.0) * effectiveTemperature * convection_factor;

        if (rise < 0.0) rise = 0.0;

        return rise;
    }

    /**
     * The surface temperature passed in is in units of Kelvin. The cloud adjustment is the
     * fraction of cloud cover obscuring each of the three major components of albedo that lie
     * below the clouds.
     * @return
     */
    public double planetAlbedo() {
        double rockFraction, cloudAdjustment, components = 0.0, cloudPart, rockPart, waterPart, icePart;

        rockFraction = 1.0 - this.hydrosphere - this.iceCover;
        if(this.hydrosphere > 0.0) {
            components += 1.0;
        }
        if(this.iceCover > 0.0) {
            components += 1.0;
        }
        if(rockFraction > 0.0) {
            components += 1.0;
        }

        cloudAdjustment = this.cloudCover / components;

        if(rockFraction >= cloudAdjustment) {
            rockFraction -= cloudAdjustment;
        } else {
            rockFraction = 0.0;
        }

        if(this.hydrosphere > cloudAdjustment) {
            this.hydrosphere -= cloudAdjustment;
        } else {
            this.hydrosphere = 0.0;
        }

        if(this.iceCover > cloudAdjustment) {
            this.iceCover -= cloudAdjustment;
        } else {
            this.iceCover = 0.0;
        }

        cloudPart = this.cloudCover * CLOUD_ALBEDO;

        if(this.surfacePressure == 0.0) {
            rockPart = rockFraction * ROCKY_AIRLESS_ALBEDO;
            icePart = this.iceCover * AIRLESS_ICE_ALBEDO;
            waterPart = 0.0;
        } else {
            rockPart = rockFraction * ROCKY_ALBEDO;
            waterPart = this.hydrosphere * WATER_ALBEDO;
            icePart = this.iceCover * ICE_ALBEDO;
        }

        return cloudPart + rockPart + waterPart + icePart;
    }

    /**
     * This function returns the dimensionless quantity of optical depth, which is useful in determining the amount of greenhouse effect on a planet.
     * @return
     */
    public double opticalDepth() {
        double opticalDepth;

        opticalDepth = 0.0;
        if(this.minimumMolecularWeight >= 0.0 && this.minimumMolecularWeight < 10.0) {
            opticalDepth += 3.0;
        }
        if(this.minimumMolecularWeight >= 10.0 && this.minimumMolecularWeight < 20.0) {
            opticalDepth += 2.34;
        }
        if(this.minimumMolecularWeight >= 20.0 && this.minimumMolecularWeight < 30.0) {
            opticalDepth += 1.0;
        }
        if(this.minimumMolecularWeight >= 30.0 && this.minimumMolecularWeight < 45.0) {
            opticalDepth += 0.15;
        }
        if(this.minimumMolecularWeight >= 45.0 && this.minimumMolecularWeight < 100.0) {
            opticalDepth += 0.05;
        }

        if(this.surfacePressure >= 70.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
            opticalDepth = opticalDepth * 8.333;
        } else if(this.surfacePressure >= 50.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
            opticalDepth = opticalDepth * 6.666;
        } else if(this.surfacePressure >= 30.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
            opticalDepth = opticalDepth * 3.333;
        } else if(this.surfacePressure >= 10.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
            opticalDepth = opticalDepth * 2.0;
        } else if(this.surfacePressure >= 5.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
            opticalDepth = opticalDepth * 1.5;
        }

        return opticalDepth;
    }

    /**
     * calculates the number of years it takes for 1/e of a gas to escape from a planet's
     * atmosphere. Taken from Dole p. 34. He cites Jeans (1916) & Jones (1923)
     * @param molecularWeight The molecular weight of the gas in question
     * @return
     */
    public double gasLife(double molecularWeight) {
        double v = rmsVelocity(molecularWeight);
        double g = this.surfaceGravity * EARTH_ACCELERATION;
        double r = this.radius * CM_PER_KM;
        double t = (Math.pow(v, 3.0) / (2.0 * Math.pow(g, 2.0) * r)) * Math.exp((3.0 * g * r) / Math.pow(v, 2.0));
        return t / (3600.0 * 24.0 * 365.256); // convert seconds to years
    }

    public double minimumMolecularWeight(Star primary) {
        double
            target = primary.age,
            guess1 = molecularLimit(),
            guess2 = guess1,
            life = gasLife(guess1);

        int	loops = 0;
        if(life > target) {
            while(life > target && loops++ < 25) {
                guess1 = guess1 / 2.0;
                life = gasLife(guess1);
            }
        } else {
            while(life < target && loops++ < 25) {
                guess2 = guess2 * 2.0;
                life = gasLife(guess2);
            }
        }

        loops = 0;
        while(guess2 - guess1 > 0.1 && loops++ < 25) {
            double guess3 = (guess1 + guess2) / 2.0;
            life = gasLife(guess3);

            if(life < target) {
                guess1 = guess3;
            } else {
                guess2 = guess3;
            }
        }

        return guess2;
    }


    /*--------------------------------------------------------------------------*/
    /*	 The temperature calculated is in degrees Kelvin.						*/
    /*	 Quantities already known which are used in these calculations:			*/
    /*		 planet->molec_weight												*/
    /*		 planet->surf_pressure												*/
    /*		 R_ecosphere														*/
    /*		 planet->a															*/
    /*		 planet->volatile_gas_inventory										*/
    /*		 planet->radius														*/
    /*		 planet->boil_point													*/
    /*--------------------------------------------------------------------------*/

    public void calculateSurfaceTemperature(Star primary, boolean first, double lastWater, double lastClouds, double lastIce, double lastTemperature, double lastAlbedo) {
        double effectiveTemperature;
        double waterRaw; // this is used for logging purposes.
        double cloudsRaw; // this is used for logging purposes.
        double greenhouseTemperature;
        boolean boilOff = false;

        if(first) {
            this.albedo = EARTH_ALBEDO;
            effectiveTemperature = effectiveTemperature(primary, this.albedo);
            greenhouseTemperature = greenhouseRise(effectiveTemperature);
            this.surfaceTemperature = effectiveTemperature + greenhouseTemperature;
            setTemperatureRange();
        }

        if(this.greenhouseEffect && this.maxTemperature < this.boilingPoint) {
            /* TODO: Logger singleton bruh
            if (flag_verbose & 0x0010)
                fprintf (stderr, "Deluge: %s %d max (%Lf) < boil (%Lf)\n",
                    planet->sun->name,
                    planet->planet_no,
                    planet->max_temp,
                    planet->boil_point);
            */
            this.greenhouseEffect = false;
            this.volatileGasInventory = volatileGasInventory(primary);
            this.surfacePressure = pressure();
            this.boilingPoint = boilingPoint();
        }

        waterRaw = this.hydrosphere = waterCoverage();
        cloudsRaw = this.cloudCover = cloudFraction();
        this.iceCover = iceFraction();

        if(this.greenhouseEffect && this.surfacePressure > 0.0) {
            this.cloudCover = 1.0;
        }

        if(this.highTemperature >= this.boilingPoint && !first && !(secondsToHoursRounded(this.dayLength) == secondsToHoursRounded(this.orbitalPeriod) || this.resonantPeriod)) {
            this.hydrosphere	= 0.0;
            boilOff = true;

            if(this.minimumMolecularWeight > WATER_VAPOR) {
                this.cloudCover = 0.0;
            } else {
                this.cloudCover = 1.0;
            }
        }

        if(this.surfaceTemperature < FREEZING_POINT_OF_WATER - 3.0) {
            this.hydrosphere = 0.0;
        }

        this.albedo = planetAlbedo();
        effectiveTemperature = effectiveTemperature(primary, this.albedo);
        greenhouseTemperature = greenhouseRise(effectiveTemperature);
        this.surfaceTemperature = effectiveTemperature + greenhouseTemperature;

        if(!first) {
            if(!boilOff) {
                this.hydrosphere = (this.hydrosphere + (lastWater * 2)) / 3;
            }
            this.cloudCover = (this.cloudCover + (lastClouds * 2)) / 3;
            this.iceCover = (this.iceCover + (lastIce * 2)) / 3;
            this.albedo = (this.albedo + (lastAlbedo * 2)) / 3;
            this.surfaceTemperature = (this.surfaceTemperature + (lastTemperature * 2)) / 3;
        }

        setTemperatureRange();

        /* TODO: Logger singleton
        if (flag_verbose & 0x0020)
            fprintf (stderr, "%5.1Lf AU: %5.1Lf = %5.1Lf ef + %5.1Lf gh%c "
                "(W: %4.2Lf (%4.2Lf) C: %4.2Lf (%4.2Lf) I: %4.2Lf A: (%4.2Lf))\n",
                planet->a,
                planet->surf_temp - FREEZING_POINT_OF_WATER,
                effective_temp - FREEZING_POINT_OF_WATER,
                greenhouse_temp,
                (planet->greenhouse_effect) ? '*' :' ',
                planet->hydrosphere, water_raw,
                planet->cloud_cover, clouds_raw,
                planet->ice_cover,
                planet->albedo);
        */
    }

    public void iterateSurfaceTemperature(Star primary) {
        double initialTemperature = estimatedTemperature(primary);

        /* TODO: This stuff might be better added to a logger singleton.
        double h2Life = gasLife(MOLECULAR_HYDROGEN);
        double h2oLife = gasLife(WATER_VAPOR);
        double n2Life = gasLife(MOLECULAR_NITROGEN);
        double nLife = gasLife(ATOMIC_NITROGEN);

        if(flag_verbose & 0x20000)
            fprintf (stderr, "%d:                     %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
                planet->planet_no,
                initial_temp,
                planet->sun->r_ecosphere, planet->a, planet->albedo
            );

        if (flag_verbose & 0x0040)
            fprintf (stderr, "\nGas lifetimes: H2 - %Lf, H2O - %Lf, N - %Lf, N2 - %Lf\n",
                h2_life, h2o_life, n_life, n2_life);
        */

        calculateSurfaceTemperature(primary, true, 0.0, 0.0, 0.0, 0.0, 0.0);

        for(int count = 0; count <= 25; count++) {
            double lastWater = this.hydrosphere;
            double lastClouds = this.cloudCover;
            double lastIce = this.iceCover;
            double lastTemperature = this.surfaceTemperature;
            double lastAlbedo = this.albedo;

            calculateSurfaceTemperature(primary, false, lastWater, lastClouds, lastIce, lastTemperature, lastAlbedo);

            if(Math.abs(this.surfaceTemperature - lastTemperature) < 0.25) {
                break;
            }
        }

        this.greenhouseRise = this.surfaceTemperature - initialTemperature;

        /* TODO: This stuff might be better added to a logger singleton.
        if (flag_verbose & 0x20000)
            fprintf (stderr, "%d: %5.1Lf gh = %5.1Lf (%5.1Lf C) st - %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
                planet->planet_no,
                planet->greenhs_rise,
                planet->surf_temp,
                planet->surf_temp - FREEZING_POINT_OF_WATER,
                initial_temp,
                planet->sun->r_ecosphere, planet->a, planet->albedo
            );
        */
    }

    /**
     * Inspired partial pressure, taking into account humidification of the air in the nasal
     * passage and throat This formula is on Dole's p. 14
     * @param gasPressure The ?surface pressure? of the gas in question
     * @return
     */
    public double inspiredPartialPressure(double gasPressure) {
        double pH2O = 47.0 * MMHG_TO_MILLIBARS; // Dole p. 15
        double fraction = gasPressure / this.surfacePressure;

        return (this.surfacePressure - pH2O) * fraction;
    }

    /**
     * This function uses figures on the maximum inspired partial pressures of Oxygen, other
     * atmospheric and traces gases as laid out on pages 15, 16 and 18 of Dole's Habitable Planets
     * for Man to derive breathability of the planet's atmosphere. - JLB
     * @return
     */
    protected atmosphereType breathableAtmosphere() {
        boolean oxygenOK = false;
        double
            MIN_O2_IPP = 72.0 * MMHG_TO_MILLIBARS, // Dole, p. 15
            MAX_O2_IPP = 400.0 * MMHG_TO_MILLIBARS; // Dole, p. 15

        if(this.atmosphere.size() == 0) {
            return atmosphereType.NONE;
        }

        Iterator<AtmosphericChemical> i = this.atmosphere.iterator();
        while(i.hasNext()) {
            AtmosphericChemical ac = i.next();
            double ipp = inspiredPartialPressure(ac.surfacePressure);
            if(ipp > ac.chem.maxInspiredPartialPressure) {
                return atmosphereType.POISONOUS;
            }
            if(ac.chem.atomicNumber == 8) { // Check for oxygen
                oxygenOK = (ipp >= MIN_O2_IPP && ipp <= MAX_O2_IPP);
            }
        }

        if(oxygenOK) {
            return atmosphereType.BREATHABLE;
        } else {
            return atmosphereType.UNBREATHABLE;
        }
    }

    /**
     * Functions for 'soft limiting' temperatures
     * @param x
     * @return
     */
    public double lim(double x) {
        return x / Math.sqrt(Math.sqrt(1.0 + x * x * x * x));
    }

    /**
     * Functions for 'soft limiting' temperatures
     * @param v
     * @param max
     * @param min
     * @return
     */
    public double soft(double v, double max, double min) {
        double dv = v - min;
        double dm = max - min;
        return (lim(2.0 * dv / dm - 1.0) + 1.0) / 2.0 * dm + min;
    }

    public void setTemperatureRange() {
        double pressmod = 1.0 / Math.sqrt(1.0 + 20.0 * this.surfacePressure / 1000.0);
        double ppmod = 1.0 / Math.sqrt(10.0 + 5.0 * this.surfacePressure / 1000.0);
        double tiltmod = Math.abs(Math.cos(this.axialTilt * Math.PI / 180.0) * Math.pow(1.0 + this.eccentricity, 2.0));
        double daymod   = 1.0 / (200.0 / this.dayLength + 1);
        double mh = Math.pow(1.0 + daymod, pressmod);
        double ml = Math.pow(1.0 - daymod, pressmod);
        double hi = mh * this.surfaceTemperature;
        double lo = ml * this.surfaceTemperature;
        double sh = hi + Math.pow((100.0 + hi) * tiltmod, Math.sqrt(ppmod));
        double wl = lo - Math.pow((150.0 + lo) * tiltmod, Math.sqrt(ppmod));
        double max = this.surfaceTemperature + Math.sqrt(this.surfaceTemperature) * 10.0;
        double min = this.surfaceTemperature / Math.sqrt(this.dayLength + 24.0);

        if(lo < min) {
            lo = min;
        }
        if(wl < 0.0) {
            wl = 0.0;
        }

        this.highTemperature = soft(hi, max, min);
        this.lowTemperature  = soft(lo, max, min);
        this.maxTemperature  = soft(sh, max, min);
        this.minTemperature = soft(wl, max, min);
    }

    protected class AtmosphericChemical {
        protected Chemical chem = null;
        protected double surfacePressure = 0.0;

        public AtmosphericChemical(Chemical c, double s) {
            this.chem = c;
            this.surfacePressure = s;
        }
    }

    public String planetType() {
        switch(this.type) {
            case tSubSubGasGiant:
                return "semi Gas Giant";
            case tSubGasGiant:
                return "sub Gas Giant";
            case tGasGiant:
                return "Gas Giant";
            case tRock:
                return "Rock";
            case tVenusian:
                return "Venusian";
            case tTerrestrial:
                return "Terrestrial";
            case tMartian:
                return "Martian";
            case tWater:
                return "Water";
            case tIce:
                return "Ice";
            case tAsteroids:
                return "Asteroids";
            case tTidallyLocked:
                return "Tidally Locked";
            default:
                return "Unknown";
        }
    }

    public String atmosphereType() {
        switch(breathableAtmosphere()) {
            case BREATHABLE:
                return "Breathable";
            case UNBREATHABLE:
                return "Unbreathable";
            case POISONOUS:
                return "Poisonous";
            default:
                return "None";
        }
    }

    public String toString() {
        if(this.gasGiant) {
            return planetType() + " " + String.format("%1$,.2f", this.sma) + "AU (" + String.format("%1$,.2f", massInJupiterMasses()) + " jm, " + String.format("%1$,.2f", massInEarthMasses()) + " em) " + numberOfMoons() + " moons" + (this.habitableMoon ? " (habitable)" : "");
        } else {
            if(this.habitable) {
                if(this.earthlike) {
                    return planetType() + " (earthlike) " + String.format("%1$,.2f", this.sma) + "AU (" + String.format("%1$,.2f", massInEarthMasses()) + " em) " + numberOfMoons() + " moons" + (this.habitableMoon ? " (habitable)" : "");
                } else {
                    return planetType() + " (habitable) " + String.format("%1$,.2f", this.sma) + "AU (" + String.format("%1$,.2f", massInEarthMasses()) + " em) " + numberOfMoons() + " moons" + (this.habitableMoon ? " (habitable)" : "");
                }
            } else {
                return planetType() + " (" + atmosphereType() + ") " + String.format("%1$,.2f", this.sma) + "AU (" + String.format("%1$,.2f", massInEarthMasses()) + " em) " + numberOfMoons() + " moons" + (this.habitableMoon ? " (habitable)" : "");
            }
        }
    }
}
