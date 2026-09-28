from ..utils.distance import euclidean
from collections import defaultdict
from .nematicPy import get_defects
from pathlib import Path
import json

class DefectTracker:
    
    def __init__(self, max_dist):
        self.defects = []
        self.active = defaultdict(list)
        self.future = defaultdict(list)
        self.max_dist = max_dist
    
    def __call__(self, ar, save=False):
        self.LX, self.LY, self.ninfo = ar.LX, ar.LY, ar.ninfo
        for frame in ar.read_frames(): self.add_frame(frame)
        if save: self.save(Path(ar._name).parent)
        else: return self.defects
    
    def add_frame(self, frame):
        for defect in get_defects(frame.QQxx, frame.QQyx, frame.LX, frame.LY):
            self.link_defect(defect, frame.time)
        # save list of acive defects for next run (all the remaining defects are inactive)
        self.active = self.future.copy()
        self.future.clear()
    
    def link_defect(self, defect, time):
        idx_dist = self.fetch_idx(defect)
        if idx_dist:# (i.e. idx is not None)
            idx, dist = idx_dist
            # link defect to stored defect
            self.defects[idx]['pos'].append(defect['pos'])
            self.defects[idx]['time'] += 1
            self.defects[idx]['distance'] += dist
            # register as active for next round
            self.future[defect['charge']].append(idx)
            # delete from current list (can't be linked more than once)
            self.active[defect['charge']].remove(idx)
        else:# we have found a new defect
            self.defects.append({
                'charge': defect['charge'],
                'pos': [defect['pos']],
                'time': 0,
                'distance': 0
            })
            # register as active for next round
            self.future[defect['charge']].append(len(self.defects) - 1)
    
    def fetch_idx(self, defect):
        """ returns idx and distance """
        distances = self.get_distances(defect)
        if distances:   # if not empty
            v, i = min([(v, i) for i, v in enumerate(distances)])
            if v < self.max_dist:
                return self.active[defect['charge']][i], v
        return None
    
    def get_distances(self, d):
        """ distance to all active defects w/ same charge """
        return [
            euclidean(
                d['pos'], self.defects[idx]['pos'][-1], [self.LX, self.LY]
            ) for idx in self.active.get(d['charge'], [])
        ]
    
    def save(self, outpath):
        with open(Path.joinpath(outpath, 'defect_trajectories.json'), 'w') as f_json:
            json.dump(self.defects, f_json)


def track_defects(ar, max_dist=5):
    defect_tracker = DefectTracker(max_dist)
    return defect_tracker(ar, save=False)

def track_and_save(ar, max_dist=5):
    defect_tracker = DefectTracker(max_dist)
    defect_tracker(ar, save=True)

