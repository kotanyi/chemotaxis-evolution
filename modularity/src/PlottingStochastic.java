import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;

public class PlottingStochastic extends JPanel {
    protected Ellipse2D.Double bug;
    protected GeneralPath gp;
    private List<Line2D.Double> scale;
    private double oldMin;
    protected double oldX, bugSize;
    protected int[] repaintBug, repaintFood;

    public PlottingStochastic(double bugPosition, List<Double> foodDistribution) {
        this.setBackground(Color.WHITE);
        bugSize = 10;
        bug = new Ellipse2D.Double(bugPosition * (Parameters.plottingWidth / Parameters.worldLength), Parameters.plottingHeight - 2 * bugSize, bugSize, bugSize);
        gp = new GeneralPath(GeneralPath.WIND_EVEN_ODD, 1000);
        gp.moveTo(0, Parameters.plottingHeight - 1.5 * bugSize - foodDistribution.get(0) * ((Parameters.plottingHeight - 3 * bugSize) / Parameters.plottingHowMuchFoodIsVisible));

        for(int i = 0; i < foodDistribution.size(); i++) {
            gp.lineTo((i / Parameters.foodEstimatesPerUnitWorldLength) * (Parameters.plottingWidth / Parameters.worldLength), Parameters.plottingHeight - 1.5 * bugSize - foodDistribution.get(i) * ((Parameters.plottingHeight - 3 * bugSize) / Parameters.plottingHowMuchFoodIsVisible));
        }

        scale = new ArrayList<Line2D.Double>();

        for(double d = Parameters.plottingHeight - 1.5 * bugSize; d >= 1.5 * bugSize; d -= (Parameters.plottingHeight - 3 * bugSize) / Parameters.plottingHowMuchFoodIsVisible) {
            scale.add(new Line2D.Double(0, d, Parameters.plottingWidth, d));
        }

        repaintBug = new int[4];
        repaintFood = new int[4];
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D graphics2D = (Graphics2D) g;

        for(int i = 0; i < scale.size(); i++) {
            if(i % 10 == 0) continue;
            graphics2D.draw(scale.get(i));
            if(i % 10 == 0) graphics2D.setColor(Color.BLACK);
        }

        graphics2D.draw(gp);

        if(bug.x == oldX) {
            graphics2D.setColor(Color.RED);
        }

        graphics2D.draw(bug);
        graphics2D.fill(bug);
        graphics2D.setColor(Color.BLACK);
    }

    public void moveBug(double bugPosition) {
        oldX = bug.x;
        bug.x = bugPosition * (Parameters.plottingWidth / Parameters.worldLength);
        repaintBug[0] = bug.x >= oldX ? (int) (oldX - 3) : (int) (bug.x - 3);
        repaintBug[1] = (int) (bug.y - 3);
        repaintBug[2] = (int) (Math.abs(bug.x - oldX) + bugSize + 6);
        repaintBug[3] = (int) (bugSize + 6);
    }

    public void moveBugRepaint() {
        repaint(repaintBug[0], repaintBug[1], repaintBug[2], repaintBug[3]);
    }

    public void moveFood(List<Double> foodDistribution) {
        gp.reset();
        gp.moveTo(0, Parameters.plottingHeight - 1.5 * bugSize - foodDistribution.get(0) * ((Parameters.plottingHeight - 3 * bugSize) / Parameters.plottingHowMuchFoodIsVisible));
        double min = Double.MAX_VALUE, y;


        for(int i = 0; i < foodDistribution.size(); i++) {
            y = Parameters.plottingHeight - 1.5 * bugSize - foodDistribution.get(i) * ((Parameters.plottingHeight - 3 * bugSize) / Parameters.plottingHowMuchFoodIsVisible);
            gp.lineTo((i / Parameters.foodEstimatesPerUnitWorldLength) * (Parameters.plottingWidth / Parameters.worldLength), y);

            if(y < min) {
                min = y;
            }
        }

        repaintFood[0] = 0;
        repaintFood[1] = (int) (min <= oldMin ? min : oldMin) - 1;
        repaintFood[2] = (int) Parameters.plottingWidth;
        repaintFood[3] = (int) (Parameters.plottingHeight - 1.5 * bugSize - repaintFood[1] + 1);
        oldMin = min;
    }

    public void moveFoodRepaint() {
        repaint(repaintFood[0], repaintFood[1], repaintFood[2], repaintFood[3]);
    }

    public Dimension getPreferredSize() {
        return new Dimension((int) Parameters.plottingWidth, (int) Parameters.plottingHeight);
    }
}
