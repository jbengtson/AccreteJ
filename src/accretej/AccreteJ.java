package accretej;

public class AccreteJ {
    public static void main(String[] args) {
        boolean habitable = false;

        /*
        System s = new System(true, false, false);
        java.lang.System.out.println(s);
        */

        System s;
        int count = 0;
        while(!habitable) {
            s = new System(true, false, false);
            if(s.habitable) {
                habitable = true;
                java.lang.System.out.println("Discarded " + count + " systems finding this one.");
                java.lang.System.out.println(s);
            } else {
                count++;
            }
        }
    }
}
