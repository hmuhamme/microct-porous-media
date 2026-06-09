from microctpm.pore_metrics import load_and_filter_image, compute_metrics

image = load_and_filter_image("sample.tif")
metrics = compute_metrics(image)

print("Pore metrics completed.")
print(metrics)
