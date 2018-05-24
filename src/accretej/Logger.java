package accretej;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Logger {
    private static Logger instance;
    private static BufferedWriter writer;
    private static SimpleDateFormat format;
    private static String LOGFILE = "accretej.log";

    public enum level {
        systemOut,
        verbose,
        specific,
        error
    }

    private Logger() {
        try {
            writer = new BufferedWriter(new FileWriter(LOGFILE, true));
            format = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");

            // initialize the log. Since we're appending, make a clear break;
            writer.newLine(); writer.newLine(); writer.newLine();
            writer.write("###############################################################################");
            writer.newLine();
            writer.write("AccreteJ log opened on " + timestamp());
            writer.newLine(); writer.newLine(); writer.newLine();
        } catch(IOException e) {
            java.lang.System.out.println("There was a non-fatal issue opening the log file: " + LOGFILE);
            e.printStackTrace();
        }
    }

    public static synchronized Logger instance() {
        if(instance == null) {
            instance = new Logger();
        }
        return instance;
    }

    public void log(level lvl, String msg) {
        String prepend = "[" + timestamp() + "]";
        try {
            switch(lvl) {
                case systemOut:
                    java.lang.System.out.println(prepend + " (OUT) - "  + msg);
                    break;
                case verbose:
                    // Verbose messages go to system out as well as the log file.
                    java.lang.System.out.println(prepend + " (VERBOSE) - "  + msg);
                    writer.write(prepend + " (VERBOSE) - "  + msg);
                    writer.newLine();
                    break;
                case specific:
                    writer.write(prepend + " (SPECIFIC) - "  + msg);
                    writer.newLine();
                    break;
                case error:
                    // Errors should go to the system error out as well as the log
                    java.lang.System.err.println(prepend + " (ERROR) - "  + msg);
                    writer.write(prepend + " (ERROR) - "  + msg);
                    writer.newLine();
                    break;
                default:
                    // report a bad enum?
                    java.lang.System.out.println("Bad enum supplied for log level.");
                    break;
            }
            writer.write(msg);
        } catch(IOException e) {
            java.lang.System.out.println("Could not write to the log file! " + LOGFILE);
            e.printStackTrace();
        }
    }

    public void close() {
        try {
            writer.flush();
            writer.close();
        } catch(IOException e) {
            java.lang.System.out.println("Could not close the log file! " + LOGFILE);
            e.printStackTrace();
        }
    }

    public String timestamp() {
        return format.format(Calendar.getInstance().getTime());
    }
}
