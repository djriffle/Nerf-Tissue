from matplotlib.cbook import flatten
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#machine learning
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split

def load_data():
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    return adata


def sc_computation(adata):
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

def gene_count_visualization(adata):
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # #plot our data
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
    sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])


def umap_and_spatial_visaulziation():

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    sc.tl.leiden(adata, key_added="clusters")

    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

    plt.show()

def positionalencoding2d(d_model, height, width):
    """
    :param d_model: dimension of the model
    :param height: height of the positions
    :param width: width of the positions
    :return: d_model*height*width position matrix
    """
    if d_model % 4 != 0:
        raise ValueError("Cannot use sin/cos positional encoding with "
                         "odd dimension (got dim={:d})".format(d_model))
    pe = np.zeros((d_model, height, width))
    # Each dimension use half of d_model
    d_model = int(d_model / 2)
    div_term = np.exp(np.arange(0., d_model, 2) *
                         -(np.log(10000.0) / d_model))
    pos_w = np.expand_dims(np.arange(0., width), axis=1)
    pos_h = np.expand_dims(np.arange(0., height), axis=1)
    pe[0:d_model:2, :, :] = np.tile(np.expand_dims(np.swapaxes(np.sin(pos_w * div_term), 0, 1), axis=1), (1, height, 1))
    pe[1:d_model:2, :, :] = np.tile(np.expand_dims(np.swapaxes(np.cos(pos_w * div_term), 0, 1), axis=1), (1, height, 1))
    pe[d_model::2, :, :] = np.tile(np.expand_dims(np.swapaxes(np.sin(pos_h * div_term), 0, 1), axis=2), (1, 1, width))
    pe[d_model + 1::2, :, :] = np.tile(np.expand_dims(np.swapaxes(np.cos(pos_h * div_term), 0, 1), axis=2), (1, 1, width))

    return pe

def machine_learning(adata):
    '''
    returns 2d positonal encodings just so we have them 
    '''
    # Turn coordinates into positional encodings
    coordinates = (adata.obsm["spatial"] / 100).astype(int)

    print("COORDINATES")
    print(coordinates)

    max_range = coordinates.max(axis=0) #makes it so we aren't dealing with offsets/blank tissue points 
    min_range = coordinates.min(axis=0) #makes it so we aren't dealing with offsets/blank tissue points

    print(f"MAX_RANGE: {max_range}")
    print(f"MIN_RANGE: {min_range}")

    #--Look into this deeper
    encoding_dim = 128
    p_enc_2d = positionalencoding2d(encoding_dim, int(max_range[0] - min_range[0]) + 1, int(max_range[1] - min_range[1]) + 1)

    print("Pos Encoding shape")
    print(p_enc_2d.shape)
    
    X = np.zeros((coordinates.shape[0], encoding_dim))
    for idx, i in enumerate(coordinates):
        X[idx] = p_enc_2d[:, i[0] - min_range[0], i[1] - min_range[1]] #puts everything back at 0 so we start at 0 index
    
    print(X[0])
    print("X shape")
    print(X.shape)
    #---

    # Gene expression
    y = adata[:, adata.var.highly_variable].X.A

    print("y shape")
    print(y.shape)

    # Split testing and training
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.05, random_state=42)

    model = keras.Sequential(
        [
            keras.Input(shape=(encoding_dim,)),
            layers.Dense(128, activation="relu"),
            layers.Dense(128, activation="relu"),
            layers.Dense(128, activation="relu"),
            layers.Dense(y.shape[1])
        ]
    )
    model.summary()

    # Call model on a test input
    opt = keras.optimizers.Adam(learning_rate=0.001)
    model.compile(optimizer=opt, loss="mse")
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)
    model.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=100, batch_size=8, callbacks=[callback])

    print("training done")

    x_cont = np.arange(0, p_enc_2d.shape[1]) + min_range[0]
    y_cont = np.arange(0, p_enc_2d.shape[2]) + min_range[0]
    xx, yy = np.meshgrid(x_cont, y_cont)
    adata

    flatten_encoding = np.transpose(p_enc_2d.reshape(encoding_dim, -1))
    print("fallten encoding shape")
    print(flatten_encoding.shape)
    print("fallten enconding")
    print(flatten_encoding)
    output = model.predict(flatten_encoding)
    print("output shape")
    print(output.shape)
    Z = output.reshape(xx.shape[1], xx.shape[0], output.shape[1]) #un  flattens output into graphable tissue
    Z = np.flip(np.swapaxes(Z, 0, 1), 0)
    print()
    print()
    print("Z shape")
    print(Z.shape)

    genes = adata.var.index[adata.var.highly_variable]
    gene_name = genes[3]
    sc.pl.spatial(adata, img_key="hires", color=[gene_name])

    plt.figure()
    plt.contour(xx, yy, Z[:, :, np.where(genes == gene_name)[0][0]], 20, cmap='jet');
    plt.show()

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(xx, yy, 
                        Z[:, :, np.where(genes == gene_name)[0][0]], 
                        cmap='jet',
                        linewidth=0, 
                        antialiased=False)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(gene_name)



    plt.show()

    # Jacobians
    x = tf.convert_to_tensor(flatten_encoding)
    with tf.GradientTape() as g:
        g.watch(x)
        with tf.GradientTape() as gg:
            gg.watch(x)
            y = model(x)
        dy_dx = gg.gradient(y, x)  # dy_dx = 2 * x
    d2y_dx2 = g.gradient(dy_dx, x)  # d2y_dx2 = 2
    print(dy_dx.shape)
    print(d2y_dx2.shape)

if __name__ == "__main__":
    print("starting")

    #sc.logging.print_versions()
    sc.set_figure_params(facecolor="white", figsize=(8, 8))
    sc.settings.verbosity = 0 # three will show hint 4 shows debug

    adata = load_data()

    sc_computation(adata)

    machine_learning(adata)


   


