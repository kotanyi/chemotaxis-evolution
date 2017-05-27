import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;

public class PlottingSingleBug extends JPanel {
    private Ellipse2D.Double bug;
    private int bugSize = 10;
    private List<Gaussian> gaussians;
    private List<GeneralPath> gaussiansPaths;

    public PlottingSingleBug(double bugPosition, List<Gaussian> gaussians) {
        this.setBackground(Color.WHITE);
        bug = new Ellipse2D.Double(bugPosition * (Parameters.plottingWidth / Parameters.worldLength), Parameters.plottingHeight - 2 * bugSize, bugSize, bugSize);
        this.gaussians = gaussians;
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D graphics2D = (Graphics2D) g;
        graphics2D.draw(bug);
        graphics2D.setColor(Color.BLACK);
        graphics2D.fill(bug);
        gaussiansPaths = new ArrayList<GeneralPath>();
        GeneralPath gp;
        double heightOfHighestGaussian = 0;

        for(Gaussian h : gaussians) {
            if(h.getHeight() > heightOfHighestGaussian) {
                heightOfHighestGaussian = h.getHeight();
            }
        }

        for(Gaussian h : gaussians) {
            gp = new GeneralPath(GeneralPath.WIND_EVEN_ODD, 1000);
            gp.moveTo(0, Parameters.plottingHeight - 1.5 * (double) bugSize);

            for(double d = 0; d <= Parameters.worldLength; d += Parameters.worldLength/1000) {
                gp.lineTo(d * (Parameters.plottingWidth / Parameters.worldLength), Parameters.plottingHeight - 1.5 * (double) bugSize - h.getFoodAt(d) * ((Parameters.plottingHeight - 3 * bugSize) / heightOfHighestGaussian));
            }

            graphics2D.draw(gp);
            gaussiansPaths.add(gp);
        }
    }

    public void moveBug(double bugPosition) {
        bug.x = bugPosition * (Parameters.plottingWidth / Parameters.worldLength);
        repaint();
    }

    public Dimension getPreferredSize() {
        return new Dimension((int) Parameters.plottingWidth, (int) Parameters.plottingHeight);
    }
}
