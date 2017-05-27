import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.font.LineMetrics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.text.DecimalFormat;
import java.text.NumberFormat;

public class SaveImages extends PlottingStochastic {
    private BufferedImage b;
    private Graphics2D g;
    private int fileCounter, textHeight, plotLabelX, timeLengthX, Y;
    private String timeLength, counter;
    FontMetrics fm;

    public SaveImages(double bugPosition, List<Double> foodDistribution, double t) {
        super(bugPosition, foodDistribution);
        fileCounter = 1;
        b = new BufferedImage((int) Parameters.plottingWidth, (int) Parameters.plottingHeight, BufferedImage.TYPE_INT_RGB);
        g = b.createGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, (int) Parameters.plottingWidth, (int) Parameters.plottingHeight);
        Font f = new Font("Dialog", Font.PLAIN, 26);
        g.setFont(f);
        fm = g.getFontMetrics();

        NumberFormat nf = new DecimalFormat("0");
        String time = nf.format(1 / Parameters.stochasticRate);
        String length = nf.format(Parameters.worldLength / Parameters.numberOfModes);
        timeLength = "Correlation time = " + time + ", correlation length = " + length;
        timeLengthX = (int) ((Parameters.plottingWidth - fm.stringWidth(timeLength)) / 2);

        LineMetrics lm = f.getLineMetrics(timeLength, g.getFontRenderContext());
        textHeight = (int) lm.getAscent();
        Y = (int) (1.5 * bugSize + textHeight);

        counter = "t = " + (int) t;

        plotLabelX = (int) Parameters.plottingWidth - textHeight - fm.stringWidth(Parameters.plotLabel);
    }

    public void changeTime(double t) {
        counter = "t = " + (int) t;
    }

    public void saveImage() {
        g.setColor(Color.WHITE);

        if(Parameters.showBug) {
            g.fillRect(repaintBug[0], repaintBug[1], repaintBug[2], repaintBug[3]);
        }

        g.fillRect(repaintFood[0], repaintFood[1], repaintFood[2], repaintFood[3]);

        if(Parameters.showCounter) {
            g.fillRect(textHeight, Y - textHeight, textHeight + fm.stringWidth(counter), Y);
        }

        g.setColor(Color.BLACK);
        g.draw(gp);

        if(Parameters.showBug) {
            g.draw(bug);
            g.fill(bug);
        }

        if(Parameters.showCorrelationTimeLength) {
            g.drawString(timeLength, timeLengthX, Y);
        }

        if(Parameters.showPlotLabel) {
            g.drawString(Parameters.plotLabel, plotLabelX, Y);
        }

        if(Parameters.showCounter) {
            g.drawString(counter, textHeight, Y);
        }

        try {
            File output = new File(System.getProperty("user.home") + "/Pictures/output" + fileCounter + ".png");
            ImageIO.write(b, "png", output);
            fileCounter++;
        } catch(IOException e) {}
    }
}
