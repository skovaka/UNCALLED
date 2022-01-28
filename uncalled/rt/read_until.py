import traceback
from threading import Event, Thread
import sys
import time
import numpy as np
import logging

try:
    import read_until
    ru_loaded = True

except ImportError:
    ru_loaded = False

if ru_loaded:
    MK_ERROR_STATUS = 0
    MK_READY = 1
    MK_STARTING = 2
    MK_PROCESSING = 3
    MK_FINISHING = 4

    RN_STARTING = 0
    RN_RUNNING = 1
    RN_FINISHING = 2
    RN_COMPLETED = 3

    class ReadUntilClient(read_until.ReadUntilClient):
        def __init__(self, 
                     mk_host=8000, 
                     mk_port="127.0.0.1", 
                     chunk_size=1.0, 
                     scan_thresh=0.99,
                     num_channels=512):

            logging.basicConfig(
                format='[%(asctime)s - %(name)s] %(message)s',
                datefmt='%H:%M:%S', level=logging.INFO)

            self.chmon_running = Event()
            self.chmon_thread = None

            read_until.ReadUntilClient.__init__(
                    self,
                    mk_host=mk_host,  
                    mk_port=mk_port,  
                    one_chunk=False,    
                    filter_strands=True,
                    prefilter_classes={'strand', 'strand1', 'strand2'}
            )


            self.in_scan = True 
            self.scan_thresh = scan_thresh
            self.chunk_size = chunk_size
            self.num_channels = num_channels

            self.ch_mux = np.zeros(num_channels, dtype=int)
            self.mux_counts = np.zeros(5, dtype=float)
            self.mux_counts[0] = len(self.ch_mux)

            self.logger.setLevel(logging.WARNING)
            self.unc_log = logging.getLogger('UNCALLED')

            #self.device = minknow.Device(self.connection)


        def run(self, steady_wait=10, scan_wait=5, refresh=0.5, **kwargs):
            self.anl_client = self.connection.analysis_configuration

            if not self._wait_for_start(steady_wait, refresh):
                return False

            self._start_chmon()

            read_until.ReadUntilClient.run(self, last_channel=self.num_channels, *kwargs)

            self.start_time = time.time()
            
            return True

        def reset(self):
            read_until.ReadUntilClient.reset(self)

            if self.chmon_thread is not None:
                self.chmon_running.clear()
                self.chmon_thread.join()
                if self.chmon_thread.is_alive():
                    self.unc_log.error("Channel monitor didn't exit right")
            self.chmon_thread = None

        def get_runtime(self):
            return time.time() - self.start_time

        def should_eject(self):
            return not self.in_scan
        
        def log(self, msg):
            self.unc_log.info(msg)

        def _start_chmon(self):
            self.start_time = time.time()
            self.last_scan = time.time()
            self.in_scan = True

            self.chmon_running.set()
            self.chmon_thread = Thread(
                target=self._run_chmon,
                name="channel_monitor"
            )
            self.chmon_thread.start()

        def _run_chmon(self, **kwargs):
            channels = self.connection.data.get_channel_states(
                first_channel = 1,
                last_channel = self.num_channels,
                use_channel_states_ids = False
            )

            try:
                self._update_muxs(channels)
            except Exception as e:
                self.unc_log.error(traceback.format_exc())

            # signal to the server that we are done with the stream.
            channels.cancel()

            self.log("Finished")

        def _scan_update(self):
            #ps = self.in_scan

            m = np.argmax(self.mux_counts)
            self.in_scan = (m != 0 and 
                            (self.mux_counts[m] / (len(self.ch_mux)-self.mux_counts[0])
                                > self.scan_thresh))

            if self.in_scan:
                self.last_scan = time.time()

            #if self.in_scan != ps:
            #    sys.stderr.write("%.2f in_scan changed to %d\n" % 
            #                     (self.get_runtime(), self.in_scan))

        def _update_muxs(self, channels):
            for channels in channels:
                if not self.chmon_running.is_set():
                    break

                for ch in channels.channel_states:
                    nmx = int(ch.config.well)
                    omx = self.ch_mux[ch.channel-1]
                    if omx != nmx:
                        self.mux_counts[omx] -= 1
                        self.mux_counts[nmx] += 1
                        self.ch_mux[ch.channel-1] = nmx

                self._scan_update()

                if self._get_minknow_status() != MK_PROCESSING:
                    self.chmon_running.clear()
                    self.running.clear()
                    break



        def _get_minknow_status(self):
            return self.connection.acquisition.current_status().status

        def _get_run_state(self):
            return self.connection.acquisition.get_acquisition_info().state


        def _wait_for_start(self, steady=10, refresh=0.01):
            if self._get_minknow_status() == MK_PROCESSING:
                self.log("Run already in progress")
                if self._update_chunk_len(False):
                    self.unc_log.error("ERROR: cannot set chunk size when run is in progress. Please stop the sequencing run and start UNCALLED before re-starting if you want to change the chunks size.")
                    return False
                return True

            processing = False
            proc_time = None
            
            self.log("Waiting for run to start")

            while True:
                status = self._get_minknow_status()
                state = self._get_run_state()

                if status == MK_STARTING or state == RN_STARTING:
                    self._update_chunk_len(True)

                if status == MK_PROCESSING:
                    if proc_time == None:
                        proc_time = time.time()

                        if not processing:
                            self.log("Waiting for steady state")
                            processing = True

                    elif time.time() - proc_time >= steady:
                        return True

                elif proc_time != None:
                    proc_time = None

                time.sleep(refresh)

            return False

        def _update_chunk_len(self, change=True):
            try:
                anl_config = self.anl_client.get_analysis_configuration()
            except:
                self.unc_log.warning("Warning: failed to check chunk size. If chunk size is set to 1 second this is fine, otherwise stop the run, use a 1 second chunk size, and please report this error.")
                return False

                if self.chunk_size != anl_config.read_detection.break_reads_after_seconds.value:
                    if change:
                        #Credit to Matt Loose for directing me to the parameter to change
                        anl_config.read_detection.break_reads_after_seconds.value=self.chunk_size
                        self.anl_client.set_analysis_configuration(anl_config)
                        self.log("Setting chunk size to %.2f" % (self.chunk_size))
                    return True
            return False


