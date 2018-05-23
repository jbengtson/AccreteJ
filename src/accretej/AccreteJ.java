package accretej;

public class AccreteJ {
    public static void main(String[] args) {
        boolean habitable = false;
        System s;
        while(!habitable) {
            s = new System(true, false, false);
            if(s.habitable) {
                habitable = true;
                java.lang.System.out.println(s);
            }
        }
    }
}
