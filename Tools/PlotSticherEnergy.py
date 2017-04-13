import sys
from PIL import Image

images1 = map(Image.open, ['42 Energy.png', '45 Energy.png', '47 Energy.png', '51 Energy.png'])
images2 = map(Image.open, ['54 Energy.png', '57 Energy.png', '60 Energy.png', '63 Energy.png'])
widths1, heights1 = zip(*(i.size for i in images1))
widths2, heights2 = zip(*(i.size for i in images2))


total_height = sum(heights1)
max_width = max(widths1)

new_im = Image.new('RGB', (max_width*2, total_height))

y_offset = 0
for im in images1:
    new_im.paste(im, (0,y_offset))
    y_offset += im.size[1]

y_offset = 0
for im in images2:
    new_im.paste(im, (max_width,y_offset))
    y_offset += im.size[1]


new_im.save('Energy.jpg')
