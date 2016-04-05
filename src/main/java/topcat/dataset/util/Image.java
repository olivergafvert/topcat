/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2016 Oliver Gäfvert

This file is part of Topcat.

Topcat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Topcat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package topcat.dataset.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.util.DistanceMatrix;
import topcat.util.SparseDistanceMatrix;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Created by oliver on 2016-03-01.
 */
public class Image {
    private static final Logger log = LoggerFactory.getLogger(Image.class);

    Pixel[][] pixels;
    public final int width;
    public final int height;

    public Image(int height, int width){
        this.width = width;
        this.height = height;
        pixels = new Pixel[height][];
        for(int i=0;i<pixels.length;i++){
            pixels[i] = new Pixel[width];
        }
    }

    public void setPixel(int x, int y, Pixel pixel){
        pixels[x][y] = pixel;
    }

    public void setPixel(int x, int y, int[] color){
        pixels[x][y] = new Pixel(color);
    }

    public Pixel getPixel(int x, int y){
        return pixels[x][y];
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(int i=0;i<height;i++){
            for(int j=0;j<width;j++){
                sb.append(Arrays.toString(pixels[i][j].color)).append(", ");
            }
            sb.append('\n');
        }
        return sb.toString();
    }

    public static DistanceMatrix getDistanceMatrix(Image image, int color, int neighbourhood){
        DistanceMatrix distanceMatrix = new SparseDistanceMatrix(image.width*image.height, image.width*image.height);

        for(int i=0;i<image.height-1;i++){
            for(int j=0;j<image.width;j++){
                int pos = i*image.width+j;
                int poscolor = image.getPixel(i, j).color[color];
                int ifrom = i-neighbourhood > 0 ? i-neighbourhood : 0;
                int ito = i + neighbourhood < image.height-1 ? i+neighbourhood : image.height-2;
                int jfrom = Math.max(j - neighbourhood, 0);
                int jto = Math.min(j + neighbourhood + 1, image.width-1);
                for(int k=ifrom;k<ito;k++){
                    for(int l=jfrom;l<jto;l++){
                        int currentColor = image.getPixel(k, l).color[color];
                        int colordiff = poscolor-currentColor;
                        colordiff = colordiff > 0 ? colordiff : -colordiff;
                        distanceMatrix.set(pos, k*image.width+l, colordiff);
                    }
                }
            }
        }
        return distanceMatrix;
    }

    private static Image loadImage(String file, Function<Integer, int[]> colourValue){
        BufferedImage img = null;
        try {
            img = ImageIO.read(new File(file));
        } catch (IOException e) {
            log.error("Failed to load image: "+file, e);
        }
        Image image = new Image(img.getHeight(), img.getWidth());
        for(int i=0;i<img.getHeight();i++){
            for(int j=0;j<img.getWidth();j++){
                image.setPixel(i, j, colourValue.apply(img.getRGB(j, i)));
            }
        }
        return image;
    }

    public static Image loadImageRGB(String file){
        return loadImage(file, Image::getRGB);
    }

    public static Image loadImageGrayscale(String file){
        return loadImage(file, Image::getGrayscale);
    }


    public static int[] getRGB(int pixel) {
        int red = (pixel >> 16) & 0xff;
        int green = (pixel >> 8) & 0xff;
        int blue = (pixel) & 0xff;
        return new int[]{red, green, blue};
    }

    public static int[] getGrayscale(int pixel) {
        int[] rgb = getRGB(pixel);
        double grayscale = 0.2989 * rgb[0] + 0.5870 * rgb[1] + 0.1140 * rgb[2];
        return new int[]{(int)Math.round(grayscale)};
    }


    public static void main(String[] args){
        Image img = Image.loadImageGrayscale("./local/image_0002.jpg");
        System.out.println(img);
//        log.info(img.width + " " + img.height);
//        double[][] d1 = Image.getDistanceMatrix(img, 0, 1);
//        double[][] d2 = Image.getDistanceMatrix(img, 1, 1);
//        for(int i=0;i<20;i++) {
//            //System.out.println(Arrays.toString(d1[i]));
//        }
//        double[] f = new double[]{0, 1};
//        PersistenceModuleCollection persistenceModuleCollection = PersistenceModuleCollection.create(new Pair<>(d1, d2),
//                new Pair<>(f, f), 2);
//        StandardNoise.computeBasicBarcode(persistenceModuleCollection.get(1).getFunctor(),
//                persistenceModuleCollection.get(1).getFiltrationValues());
//        try {
//            FFmpegFrameGrabber g = new FFmpegFrameGrabber("textures/video/anim.mp4");
//
//            g.start();
//
//            for (int i = 0; i < 50; i++) {
//                ImageIO.write(g.grabImage(), "png", new File("frame-dump/video-frame-" + System.currentTimeMillis() + ".png"));
//            }
//
//            g.stop();
//        }catch (Exception e){
//            e.printStackTrace();
//        }
    }


    public class Pixel{
        private int[] color;

        public Pixel(int[] color){
            this.color = color;
        }

        public int[] getColor(){
            return color;
        }
    }
}
