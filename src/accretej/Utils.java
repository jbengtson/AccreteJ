package accretej;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Random;
import java.util.Vector;
import java.lang.System;

class Utils {
    private static Utils instance;
    private static Random random;
    private static Vector<Star>
        MStars,
        KStars,
        GStars,
        FStars,
        AStars,
        BStars,
        OStars,
        WStars,
        CStars,
        SStars,
        White,
        Giant;
    private static Chemical[] Chemtable;

    public final static double
        ECCENTRICITY_COEFF = 0.077,
        EARTH_SURF_PRES_IN_MILLIBARS = 1013.25,
        EARTH_SURF_PRES_IN_MMHG = 760.0,
        MMHG_TO_MILLIBARS = EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_MMHG,
        PPM_PRESSURE = EARTH_SURF_PRES_IN_MILLIBARS / 1000000.0;

    private Utils() {
        // we want to initialize the random instance.
        random = new Random(java.lang.System.currentTimeMillis());
        loadChemicals();
        loadStars();
    }

    public static synchronized Utils instance() {
        if(instance == null) {
            instance = new Utils();
        }
        return instance;
    }

    public double randomNumber(double inner, double outer) {
        return random.nextDouble() * (outer - inner) + inner;
    }

    public double randomEccentricity() {
        return(1.0 - Math.pow(random.nextDouble(), ECCENTRICITY_COEFF));
    }

    public double about(double value, double variation) {
        return value + (value * randomNumber(-variation, variation));
    }

    public int countPlanets(Planet p) {
        if(p == null) {
            return 0;
        }
        Planet current = p.next;
        int x = 1;
        while(current != null) {
            x++;
            current = current.next;
        }
        return x;
    }

    public double sumMassOfPlanets(Planet p) {
        if(p == null) {
            return 0.0;
        }
        Planet current = p.next;
        double x = p.dustMass + p.gasMass;
        while(current != null) {
            x += current.dustMass + current.gasMass;
            current = current.next;
        }
        return x;
    }

    public Star randomStar() {
        double roll = random.nextDouble();
        if(roll <= 0.907) { // Main sequence stars
            roll = random.nextDouble();
            if(roll <= 0.751) { // M type main sequence stars
                return MStars.get(random.nextInt(MStars.size())).deviate();
            } else if(roll <= 0.887) { // K type main sequence stars
                return KStars.get(random.nextInt(KStars.size())).deviate();
            } else if(roll <= 0.960) { // G type main sequence stars
                return GStars.get(random.nextInt(GStars.size())).deviate();
            } else if(roll <= 0.991) { // F type main sequence stars
                return FStars.get(random.nextInt(FStars.size())).deviate();
            } else { // A type main sequence stars
                return AStars.get(random.nextInt(AStars.size())).deviate();
            }
        } else if(roll <= 0.969) { // White dwarves
            return White.get(random.nextInt(White.size())).deviate();
        } else if(roll <= 0.998) { // Giants
            return Giant.get(random.nextInt(Giant.size())).deviate();
        } else { // Other stars
            roll = random.nextDouble();
            if(roll <= 0.785) { // B type stars
                return BStars.get(random.nextInt(BStars.size())).deviate();
            } else if(roll <= 0.999) { // O type stars
                return OStars.get(random.nextInt(OStars.size())).deviate();
            } else { // Specials and "odd" stars
                roll = random.nextDouble();
                if(roll <= 0.997) { // O type stars
                    return OStars.get(random.nextInt(OStars.size())).deviate();
                } else if(roll <= 0.998) { // Wolf Rayet type stars
                    return WStars.get(random.nextInt(WStars.size())).deviate();
                } else if(roll <= 0.999) { // C type stars
                    return CStars.get(random.nextInt(CStars.size())).deviate();
                } else { // S type stars
                    return SStars.get(random.nextInt(SStars.size())).deviate();
                }
            }
        }
    }

    public Star randomMStar() {
        return MStars.get(random.nextInt(MStars.size())).deviate();
    }

    public Star randomKStar() {
        return KStars.get(random.nextInt(KStars.size())).deviate();
    }

    public Star randomGStar() {
        return GStars.get(random.nextInt(GStars.size())).deviate();
    }

    public Star randomFStar() {
        return FStars.get(random.nextInt(FStars.size())).deviate();
    }

    public Chemical[] getChemicals() {
        return Chemtable;
    }

    private static void loadStars() {
        MStars = loadStarType("data/MV_Stars.csv");
        KStars = loadStarType("data/KV_Stars.csv");
        GStars = loadStarType("data/GV_Stars.csv");
        FStars = loadStarType("data/FV_Stars.csv");
        AStars = loadStarType("data/AV_Stars.csv");
        BStars = loadStarType("data/B_Stars.csv");
        OStars = loadStarType("data/O_Stars.csv");
        WStars = loadStarType("data/WR_Stars.csv");
        CStars = loadStarType("data/C_Stars.csv");
        SStars = loadStarType("data/S_Stars.csv");
        White = loadStarType("data/White_Dwarves.csv");
        Giant = loadStarType("data/Giants.csv");
    }

    private static Vector<Star> loadStarType(String filename) {
        Vector<Star> stars = new Vector<>();
        Star s;
        String line;
        String[] split;

        try {
            BufferedReader input = new BufferedReader(new FileReader(filename));
            while((line = input.readLine()) != null) {
                split = line.split(",");
                s = new Star(Double.valueOf(split[1]), Double.valueOf(split[2]), Double.valueOf(split[3]), Double.valueOf(split[4]), Double.valueOf(split[6]));
                s.stellarType = split[0];
                s.red = Integer.valueOf(split[9]);
                s.green = Integer.valueOf(split[10]);
                s.blue = Integer.valueOf(split[11]);
                stars.add(s);
            }
            input.close();
        } catch(Exception e) {
            System.err.println(filename);
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
        return stars;
    }

    private static void loadChemicals() {
        Chemtable = new Chemical[13];
        Chemtable[0] = new Chemical(1,"H", "Hydrogen", 1.0079, 14.06, 20.40, 8.99e-05, 0.00125893, 27925.4, 1, 0.0);
        Chemtable[1] = new Chemical(2, "He", "Helium", 4.0026, 3.46, 4.20, 0.0001787, 7.94328e-09, 2722.7, 0, 61000.0 * MMHG_TO_MILLIBARS);
        Chemtable[2] = new Chemical(7, "N", "Nitrogen", 14.0067, 63.34, 77.40, 0.0012506, 1.99526e-05, 3.13329, 0, 2330.0 * MMHG_TO_MILLIBARS);
        Chemtable[3] = new Chemical(8, "O", "Oxygen", 15.9994, 54.80, 90.20, 0.001429, 0.501187, 23.8232,  10, 400.0 * MMHG_TO_MILLIBARS);
        Chemtable[4] = new Chemical(10, "Ne", "Neon", 20.1700, 24.53, 27.10, 0.0009, 5.01187e-09, 3.4435e-5, 0, 3900.0 * MMHG_TO_MILLIBARS);
        Chemtable[5] = new Chemical(18, "Ar", "Argon", 39.9480, 84.00, 87.30, 0.0017824, 3.16228e-06, 0.100925, 0, 1220.0 * MMHG_TO_MILLIBARS);
        Chemtable[6] = new Chemical(36, "Kr", "Krypton", 83.8000, 116.60, 119.70, 0.003708, 1e-10, 4.4978e-05, 0, 350.0 * MMHG_TO_MILLIBARS);
        Chemtable[7] = new Chemical(54, "Xe", "Xenon", 131.3000, 161.30, 165.00, 0.00588, 3.16228e-11, 4.69894e-06, 0, 160.0 * MMHG_TO_MILLIBARS);
        Chemtable[8] = new Chemical(900, "NH3", "Ammonia", 17.0000, 195.46, 239.66, 0.001,0.002, 0.0001, 1, 100.0 * PPM_PRESSURE);
        Chemtable[9] = new Chemical(901, "H2O", "Water",  18.0000, 273.16, 373.16, 1.000, 0.03, 0.001, 0, 0.0);
        Chemtable[10] = new Chemical(902, "CO2", "CarbonDioxide", 44.0000, 194.66, 194.66, 0.001, 0.01, 0.0005, 0, 7.0 * MMHG_TO_MILLIBARS);
        Chemtable[11] = new Chemical(903, "O3", "Ozone", 48.0000, 80.16, 161.16, 0.001, 0.001, 0.000001, 2, 0.1 * PPM_PRESSURE);
        Chemtable[12] = new Chemical(904, "CH4", "Methane", 16.0000, 90.16, 109.16, 0.010, 0.005, 0.0001, 1, 50000.0 * PPM_PRESSURE);

/* TODO: Format this correctly in case they can be added later? Some good stuff here.
        {AN_F,  "F",  "F",							 "Fluorine",        18.9984,  53.58,  85.10,  0.001696,  0.000630957, 0.000843335,   50,	MAX_F_IPP},
        {AN_CL, "Cl", "Cl",							 "Chlorine",        35.4530, 172.22, 239.20,  0.003214,  0.000125893, 0.005236,      40,	MAX_CL_IPP},
        { 910, "H2", "H2",  2, 14.06, 20.40, 8.99e-05,  0.00125893, 27925.4  },
        { 911, "N2", "N2", 28, 63.34, 77.40, 0.0012506, 1.99526e-05,3.13329  },
        { 912, "O2", "O2", 32, 54.80, 90.20, 0.001429,  0.501187, 23.8232, 10},
        {AN_CH3CH2OH,"CH3CH2OH", "Ethanol",  46.0000, 159.06, 351.66,  0.895,     0.001,       0.001,         0},
*/

    }
}
