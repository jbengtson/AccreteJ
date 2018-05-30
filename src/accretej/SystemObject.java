package accretej;

public abstract class SystemObject {
    public final static double
        G = 6.67408E-11,
        SUN_MASS = 1.989E30, // Mass of the Sun in kg
        SUN_MASS_IN_GRAMS = 1.989E33, // Mass of the Sun on g
        SUN_RADIUS = 6.96392E8, // radius of the Sun in m
        SUN_MASS_IN_EARTH_MASSES = 333054.25,
        JUPITER_MASS = 1.8982E27, // Mass of Jupiter in kg
        JUPITER_RADIUS = 6.9911E7, // radius of Jupiter in m
        EARTH_MASS = 5.97237E24, // mass of Earth in kg
        EARTH_RADIUS = 6.371E6, // radius of Earth in m
        EARTH_SMA = 149597870700.0, // for calculation of an astronomical unit
        KM_PER_AU = 149597870.7, // Both of these are set to the metric SI unit for AU
        CM_PER_KM = 100000.0,
        CM_PER_METER = 100.0,
        // MOLAR_GAS_CONST = 8314.41, // units: g*m2/(sec2*K*mol)
        MOLAR_GAS_CONST = 8.3144621, // units of J/molK
        GAS_RETENTION_THRESHOLD = 6.0,
        EARTH_ACCELERATION = 9.80655, // m/s
        EARTH_DENSITY = 5.514, // g/cc
        FREEZING_POINT_OF_WATER = 273.15, // Units of degrees Kelvin
        EARTH_AVERAGE_CELSIUS = 14.0, // Average Earth Temperature
        EARTH_AVERAGE_KELVIN = EARTH_AVERAGE_CELSIUS + FREEZING_POINT_OF_WATER,
        EARTH_WATER_MASS_PER_AREA = 3.83E15, // grams per square km
        EARTH_SURF_PRES_IN_MILLIBARS = 1013.25,
        EARTH_EXOSPHERIC_TEMPERATURE = 1273.0, // in Kelvin
        MMHG_TO_MILLIBARS = 1.33322,
        ASTEROID_MASS_LIMIT = 0.001;

    // Molecular weights seperated for because.
    public final static double
        ATOMIC_HYDROGEN = 1.0,
        MOLECULAR_HYDROGEN = 2.0,
        HELIUM = 4.0,
        ATOMIC_NITROGEN = 14.0,
        ATOMIC_OXYGEN = 16.0,
        METHANE = 16.0,
        AMMONIA = 17.0,
        WATER_VAPOR = 18.0,
        NEON = 20.2,
        MOLECULAR_NITROGEN = 28.0,
        CARBON_MONOXIDE = 28.0,
        NITRIC_OXIDE = 30.0,
        MOLECULAR_OXYGEN = 32.0,
        HYDROGEN_SULFIDE = 34.1,
        ARGON = 39.9,
        CARBON_DIOXIDE = 44.0,
        NITROUS_OXIDE = 44.0,
        NITROGEN_DIOXIDE = 46.0,
        OZONE = 48.0,
        SULFUR_DIOXIDE = 64.1,
        SULFUR_TRIOXIDE = 80.1,
        KRYPTON = 83.8,
        XENON = 131.3;

    // Albedos
    public final double
        ICE_ALBEDO = 0.7,
        CLOUD_ALBEDO = 0.52,
        GAS_GIANT_ALBEDO = 0.5,
        AIRLESS_ICE_ALBEDO = 0.5,
        EARTH_ALBEDO = 0.3,
        GREENHOUSE_TRIGGER_ALBEDO = 0.20,
        ROCKY_ALBEDO = 0.15,
        ROCKY_AIRLESS_ALBEDO = 0.07,
        WATER_ALBEDO = 0.04;

    protected double
        mass, // This is in solar masses, use a conversion method if needed.
        sma,
        asMoonSMA = 0.0,
        asMoonEccentricty = 0.0,
        eccentricity,
        inclination;

    /**
     * @return The Standard Gravitational Parameter of this object.
     */
    public double mu() {
        return massInKg() * G;
    }

    public static double massInKg(double mass) {
        return mass * SUN_MASS;
    }

    public double massInKg() {
        return massInKg(this.mass);
    }

    public static double massInGrams(double mass) {
        return mass * SUN_MASS_IN_GRAMS;
    }

    public double massInGrams() {
        return massInGrams(this.mass);
    }

    public static double massInJupiterMasses(double mass) {
        return massInKg(mass) / JUPITER_MASS;
    }

    public double massInJupiterMasses() {
        return massInJupiterMasses(this.mass);
    }

    public static double massInEarthMasses(double mass) {
        return massInKg(mass) / EARTH_MASS;
    }

    public double massInEarthMasses() {
        return massInEarthMasses(this.mass);
    }

    /**
     * @return The orbital apoapsis
     */
    public static double apoapsis(double sma, double ecc) {
        return sma * (1.0 + ecc);
    }

    /**
     * @return The orbital periapsis
     */
    public static double periapsis(double sma, double ecc) {
        return sma * (1.0 - ecc);
    }

    /**
     * @return The orbital apoapsis
     */
    public double apoapsis() {
        return apoapsis(this.sma, this.eccentricity);
    }

    /**
     * @return The orbital periapsis
     */
    public double periapsis() {
        return periapsis(this.sma, this.eccentricity);
    }

    /**
     * @return The orbital apoapsis
     */
    public double asMoonApoapsis() {
        return apoapsis(this.asMoonSMA, this.asMoonEccentricty);
    }

    /**
     * @return The orbital periapsis
     */
    public double asMoonPeriapsis() {
        return periapsis(this.asMoonSMA, this.asMoonEccentricty);
    }

    public static double smaInMeters(double sma) {
        return sma * EARTH_SMA;
    }

    public double smaInMeters() {
        return smaInMeters(this.sma);
    }

    public static double AUtoKm(double au) {
        return au * KM_PER_AU;
    }

    /**
     * Returns the period of this orbit given the gravitational parameter of the central body.
     * @param mu The Standard Gravitational Parameter of the central body. See https://en.wikipedia.org/wiki/Standard_gravitational_parameter
     * @return The orbital period in seconds.
     */
    public double orbitalPeriod(double mu) {
        return 2.0 * Math.PI * Math.sqrt(Math.pow(sma, 3.0) / mu);
    }

    public int secondsToHoursRounded(double sec) {
        return (int)(sec / 3600.0);
    }

    public double secondsToHours(double sec) {
        return sec / 3600.0;
    }

    public int secondsToDaysRounded(double sec) {
        return (int)(sec / 86400);
    }

    public double secondsToDays(double sec) {
        return sec / 86400.0;
    }

    public double secondsToYears(double sec) {
        return sec / 31557600;
    }

    /**
     * @return The orbital SMA expressed as a unit of astronomical distance (AU, one Earth SMA).
     */
    public double AU() {
        return sma / EARTH_SMA;
    }
}
