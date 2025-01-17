
# Detectar si se está ejecutando en un entorno SLURM
if (Sys.getenv("SLURM_JOB_ID") != "") {
  # Si estamos en SLURM, usar el número de núcleos especificado por SLURM
  n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK"))  # SLURM cores
} else {
  # Si estamos en una máquina local, detectar los núcleos disponibles menos 1
  n.cores <- parallel::detectCores() - 1
}

# Crear el cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

# Verificar los núcleos utilizados
print(paste("Núcleos utilizados:", n.cores))
