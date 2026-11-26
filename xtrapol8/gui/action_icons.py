"""Scalable Qt-drawn toolbar icons (no external image assets required)."""
import math
from PyQt5.QtCore import Qt, QPointF, QRectF
from PyQt5.QtGui import QColor, QIcon, QPainter, QPainterPath, QPen, QPixmap


def action_icon(action, size=64):
    pixmap = QPixmap(size, size)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    center = size / 2.0
    radius = size * 0.40
    stroke = size * 0.075
    if action == 'run':
        color = QColor('#24965A')
        painter.setPen(QPen(color, stroke))
        painter.setBrush(Qt.NoBrush)
        painter.drawEllipse(QRectF(center-radius, center-radius, 2*radius, 2*radius))
        triangle = QPainterPath()
        triangle.moveTo(size*.41, size*.31)
        triangle.lineTo(size*.70, size*.50)
        triangle.lineTo(size*.41, size*.69)
        triangle.closeSubpath()
        painter.setPen(Qt.NoPen)
        painter.setBrush(color)
        painter.drawPath(triangle)
    elif action == 'stop':
        color = QColor('#D64B4B')
        octagon = QPainterPath()
        for i in range(8):
            angle = math.radians(22.5 + 45*i)
            point = QPointF(center+radius*math.cos(angle), center+radius*math.sin(angle))
            if i == 0:
                octagon.moveTo(point)
            else:
                octagon.lineTo(point)
        octagon.closeSubpath()
        painter.setPen(QPen(color, stroke, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
        painter.setBrush(Qt.NoBrush)
        painter.drawPath(octagon)
        painter.setPen(QPen(color, stroke*1.3, Qt.SolidLine, Qt.RoundCap))
        painter.drawLine(QPointF(size*.39, size*.39), QPointF(size*.61, size*.61))
        painter.drawLine(QPointF(size*.61, size*.39), QPointF(size*.39, size*.61))
    else:
        painter.end()
        raise ValueError(f'Unknown action: {action}')
    painter.end()
    return QIcon(pixmap)
