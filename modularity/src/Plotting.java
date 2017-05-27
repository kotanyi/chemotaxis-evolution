import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.IRangePolicy;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.rangepolicies.RangePolicyFixedViewport;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import info.monitorenter.util.Range;
import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.List;

public class Plotting {
    private Chart2D chart;
    private List<ITrace2D> traces;
    private JFrame frame;

    public Plotting() {
        super();
        traces = new ArrayList<ITrace2D>();
    }

    public void initialise() {
        chart = new Chart2D();
        Range rangeY = new Range(-0.1, 1.1);
        IRangePolicy rangePolicyY = new RangePolicyFixedViewport(rangeY);
        Range rangeX = new Range(-110, 250);

        if(Parameters.stimulusCourse[0]) { rangeX = new Range(-110, 160); }
        else if(Parameters.stimulusCourse[1]) { rangeX = new Range(-110, 250); }
        else if(Parameters.stimulusCourse[2]) { rangeX = new Range(-0.1, 3.64); }
        else if(Parameters.stimulusCourse[3]) { rangeX = new Range(-110, 885); }
        else if(Parameters.stimulusCourse[4]) { rangeX = new Range(-10, 3010); }

        IRangePolicy rangePolicyX = new RangePolicyFixedViewport(rangeX);
        chart.getAxisY().setRangePolicy(rangePolicyY);
        chart.getAxisX().setRangePolicy(rangePolicyX);

        traces.add(new Trace2DSimple());
        traces.get(0).setColor(Color.cyan);
        traces.get(0).setName("Output");
        chart.addTrace(traces.get(0));

        traces.add(new Trace2DSimple());
        traces.get(1).setColor(Color.blue);
        traces.get(1).setName("Input");
        chart.addTrace(traces.get(1));

        traces.add(new Trace2DSimple());
        traces.get(2).setColor(Color.red);
        traces.get(2).setName("Sensitivity");
        //chart.addTrace(traces.get(2));

        traces.add(new Trace2DSimple());
        traces.get(3).setColor(Color.green);
        traces.get(3).setName("Goal");
        chart.addTrace(traces.get(3));

        frame = new JFrame("Modularity");
        frame.setSize(1200, 500);
        frame.getContentPane().add(chart);

        frame.addWindowListener(
            new WindowAdapter() {
                public void windowClosing(WindowEvent windowEvent) {
                    System.exit(0);
                }
            }
        );

        frame.setVisible(true);
    }

    public void clearTraces() {
        for(ITrace2D trace : traces)
            trace.removeAllPoints();
    }

    public void updateTrace(int trace, double x, double y) {
        traces.get(trace).addPoint(x, y);
    }

    public void destroy() {
        frame.dispose();
    }
}
