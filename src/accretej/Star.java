package accretej;

public class Star extends SystemObject {
    public final static double STELLAR_DEVIATION = 0.05; // Amount by which individual stars might vary from their base type

    protected String
        stellarType = "",
        name = "",
        desig = "";
    protected int
        red = 0,
        green = 0,
        blue = 0;
    protected double
        luminosity = 0.0,
        radius = 0.0,
        temperature = 0.0,
        absoluteMagnitude = 0.0,
        lifeTime = 0.0,
        age = 0.0,
        radiusEcosphere = 0.0,
        M2 = 0.0,
        A = 0.0,
        E = 0.0;

    public Star(double smass, double slum, double srad, double temp, double mag) {
        this.mass = smass;
        this.luminosity = slum;
        this.radius = srad;
        this.temperature = temp;
        this.absoluteMagnitude = mag;
        this.radiusEcosphere = Math.sqrt(luminosity); // is this accurate? Presumably this scales with the inverse square law, so it sounds right.
        this.lifeTime = 1.0E10 * (mass / luminosity);

        recalc();
    }

    public void recalc() {
        this.lifeTime = 10E10 * Math.pow(this.mass, 2.5); // http://hyperphysics.phy-astr.gsu.edu/hbase/Astro/startime.html
    }

    public void setAge() {
        if(lifeTime < 6.0E9) {
            this.age = Utils.instance().randomNumber(1.0E9, lifeTime);
        } else {
            this.age = Utils.instance().randomNumber(1.0E9, 6.0E9);
        }
    }

    public double stellarDustLimit() {
        return 200.0 * Math.pow(this.mass, 1.0 / 3.0);
    }

    public double innermostPlanet() {
        return 0.3 * Math.pow(mass, 1.0 / 3.0); // TODO: Check these numbers to ensure accuracy
    }

    public double outermostPlanet() {
        return 50.0 * Math.pow(mass, 1.0 / 3.0); // TODO: Check these numbers to ensure accuracy
    }

    /**
     * @return A copy of this star
     */
    public Star copy() {
        Star s = new Star(this.mass, this.luminosity, this.radius, this.temperature, this.absoluteMagnitude);
        s.stellarType = this.stellarType;
        s.red = this.red;
        s.green = this.green;
        s.blue = this.blue;

        return s;
    }

    /**
     * @return A copy of this star deviated by a random amount
     */
    public Star deviate() {
        Star s = this.copy();
        double v = Utils.instance().about(STELLAR_DEVIATION, 1);
        s.mass = s.mass + s.mass * v;
        s.luminosity = s.luminosity + s.luminosity * v;
        s.radius = s.radius + s.radius * v;
        s.temperature = s.temperature + s.temperature * v;
        s.recalc();

        return s;
    }

    public String toString() {
        return stellarType + " (" + String.format("%1$,.2f", mass) + "sm)";
    }
}
