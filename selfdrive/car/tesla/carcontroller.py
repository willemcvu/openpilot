from common.realtime import DT_CTRL
from common.numpy_fast import clip, interp
from selfdrive.car.tesla.teslacan import TeslaCAN
from opendbc.can.packer import CANPacker
from selfdrive.car.tesla.values import CarControllerParams

class CarController():
  def __init__(self, dbc_name, CP, VM):
    self.CP = CP
    self.last_angle = 0
    self.packer = CANPacker(dbc_name)
    self.tesla_can = TeslaCAN(dbc_name, self.packer)

    self.torque_blending_integral = 0
    self.user_torque_moving_average = 0

  def compute_torque_blending(self, driver_torque):
    self.user_torque_moving_average = ((self.user_torque_moving_average * 9) + driver_torque) / 10
    if abs(self.user_torque_moving_average) > 0.1:
      self.torque_blending_integral += driver_torque * DT_CTRL
    else:
      self.torque_blending_integral *= 0.98

    angle_offset = self.torque_blending_integral * 2 # TODO: tune this value
    return angle_offset

  def update(self, enabled, CS, frame, actuators, cruise_cancel):
    can_sends = []

    # Temp disable steering on a hands_on_fault
    hands_on_fault = (CS.steer_warning == "EAC_ERROR_HANDS_ON" and CS.hands_on_level >= 3)
    lkas_enabled = enabled and not hands_on_fault

    if lkas_enabled:
      apply_angle = actuators.steeringAngleDeg

      # Driver torque blending
      apply_angle += self.compute_torque_blending(CS.out.steeringTorque)

      # Angular rate limit based on speed
      steer_up = (self.last_angle * apply_angle > 0. and abs(apply_angle) > abs(self.last_angle))
      rate_limit = CarControllerParams.RATE_LIMIT_UP if steer_up else CarControllerParams.RATE_LIMIT_DOWN
      max_angle_diff = interp(CS.out.vEgo, rate_limit.speed_points, rate_limit.max_angle_diff_points)
      apply_angle = clip(apply_angle, (self.last_angle - max_angle_diff), (self.last_angle + max_angle_diff))
    else:
      apply_angle = CS.out.steeringAngleDeg

    self.last_angle = apply_angle
    can_sends.append(self.tesla_can.create_steering_control(apply_angle, lkas_enabled, frame))

    # Cancel when openpilot is not enabled anymore
    if not enabled and bool(CS.out.cruiseState.enabled):
      cruise_cancel = True

    if ((frame % 10) == 0 and cruise_cancel):
      can_sends.append(self.tesla_can.create_action_request(CS.msg_stw_actn_req, cruise_cancel))

    # TODO: HUD control

    return can_sends
