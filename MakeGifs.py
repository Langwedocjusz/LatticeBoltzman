import os
import imageio

def GifFromImages(img_dir, output_path):
    images = []

    names = []

    for filename in os.listdir(img_dir):
        filepath = os.path.join(img_dir, filename)

        if not os.path.isfile(filepath): 
            continue

        name, extension = filename.split('.')

        if extension == "png" and name.isnumeric():
            names.append(name)

    names = [int(name) for name in names]
    names.sort()

    for name in names:
        filepath = os.path.join(img_dir, str(name) + ".png")

        images.append(imageio.imread(filepath))

    imageio.mimsave(output_path, images, format='GIF', duration=4, loop=0)


GifFromImages('plots/density', 'plots/density.gif')
GifFromImages('plots/velx', 'plots/velx.gif')
GifFromImages('plots/vely', 'plots/vely.gif')