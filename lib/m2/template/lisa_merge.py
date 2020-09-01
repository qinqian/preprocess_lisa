import h5py
import numpy as np

with h5py.File(snakemake.input.h5[0]) as inf:
    nrp = inf["RP"].shape[0]
    nc = inf["OrderCount"].shape[0]

with h5py.File(snakemake.output.fh5, "a") as store:
    ct = store.create_dataset("OrderCount", dtype=np.float32, shape=(nc, len(snakemake.input.h5)), compression='gzip', shuffle=True, fletcher32=True)
    RP = store.create_dataset("RP", dtype=np.float32, shape=(nrp, len(snakemake.input.h5)), compression='gzip', shuffle=True, fletcher32=True)
    ids = []
    for i, d in enumerate(snakemake.input.h5):
        with h5py.File(d) as inf:
            ct[:,i] = inf["OrderCount"][:,0]
            RP[:,i] = inf["RP"][:,0]
        store.flush()
        ids.append(str.encode(snakemake.params.prefix + ".%s" % snakemake.params.labels[i], 'utf-8'))
    store["IDs"] = np.array(ids)
    store.flush()
