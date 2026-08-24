# hydrogen-model — offline JupyterLab environment for the hierarchical hydrogen permeation model.
#
# Architecture is NOT pinned here on purpose; it is a release decision:
#     docker build --platform linux/amd64 -t hydrogen-model:1.0.0-rc1 container/
#
# The application code is NOT copied in. It is bind-mounted at runtime so the
# model can be updated without rebuilding, and so one validated image serves
# successive versions of the notebooks.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Arbitrary-UID tolerance. A fixed-UID user is WORSE than root for bind mounts:
# container UID 1000 writing to a host dir owned by 1001 gets EACCES and
# notebooks cannot save at all, whereas root always writes. But `--user $(id -u)`
# alone also fails, because an unknown UID gets HOME=/ which is not writable and
# JupyterLab needs a writable HOME. So: own the home dir by group 0, make it
# group-writable, run as 1000:0, and keep Jupyter's state in /tmp. Any host UID
# then works via `--user "$(id -u):0"`.
RUN useradd --create-home --uid 1000 --gid 0 --shell /bin/bash appuser && \
    chmod -R g=u /home/appuser

# PYTHONPATH=/workspace so the bind-mounted code is importable no matter how
# python is invoked. Without it, `python /opt/model_smoke_test.py` puts /opt on
# sys.path[0] -- not the cwd -- and `import calculations` fails even with
# -w /workspace. The notebooks do their own sys.path.insert, but scripts do not.
ENV PYTHONPATH=/workspace \
    HOME=/home/appuser \
    JUPYTER_CONFIG_DIR=/tmp/.jupyter \
    JUPYTER_DATA_DIR=/tmp/.jupyter/data \
    JUPYTER_RUNTIME_DIR=/tmp/.jupyter/runtime

WORKDIR /workspace

# Install from the lock, never from requirements.in — the lock is the validated
# environment. `pip check` fails the build on a broken resolve.
COPY requirements.lock.txt /opt/requirements.lock.txt
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir -r /opt/requirements.lock.txt && \
    python -m pip check

COPY verify_environment.py model_smoke_test.py /opt/

USER 1000:0

EXPOSE 8888

# No --allow-root: this does not run as root. JupyterLab prints a tokenised URL
# on startup; the README tells the user to copy it.
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser"]
